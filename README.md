# GATKcall-with-BQSR-from-GVCF-to-VCF-（gVCF → 染色体並列 → 結合）
遺伝研スパコン上でBQSR（ベースクオリティを向上させることでsnpcallの精度を向上させる方法）を採用したGATKcallのパイプラインを掲載する.
一つ前のファイルでVCFまで作成したものを統合し、VCFファイルまで作成するまでのパイプライン。
COMBINE GVCFを用いるのではなく、GenomicsDBimportを用いて容量の軽減かつ効率化を図っている。
また、染色体分割を用いて並列処理を可能にして効率化を図っている。

このスクリプトは、ランチャー（親ジョブ）が配列ジョブ（染色体ごと）を一括投入し、すべてのワーカーが完了したら 結合ジョブ を自動起動して全染色体を連結する構成です。
再実行すると、既に出力済みの染色体VCFはスキップされ、未完了ぶんだけ進みます。

## 1. ジョブ構成（全体像）

ランチャー（このスクリプト本体）
軽量リソースで OK。ワーカー配列と結合ジョブを sbatch で投入して終了。

ワーカー（MODE=worker）
染色体ごとに GenomicsDBImport → GenotypeGVCFs を実行。並列度は --array=1-N%ARRAY_MAXPAR。

結合ジョブ（MODE=concat）
すべての染色体VCFの生成完了を前提に bcftools concat で一つに連結。

>[!NOTE]
依存関係は --dependency=afterok:<配列ジョブID> で管理。ワーカーが1つでも失敗すると結合ジョブは走りません。

## 2. 主要パラメータ（環境変数で上書き可）
```bash
# 参照・入力
genome=201123_Psic.polished.scaffold.masked.10k
REF_DIR=$HOME/working/Pipeline_SNP_call/reference
CHR_LIST=$REF_DIR/chromosomes.list     # 1行=1染色体名、空行と#は無視

# gVCF探索ルート & 実行ルート
OUT_ROOT=$HOME/working/pipeline_out     # <sample>.g.vcf.gz があること
RUN_DIR=${OUT_ROOT}/_joint_YYYYmmdd_HHMMSS  # 自動作成

# コンテナ/コマンド
GATK_SIF=/usr/local/biotools/g/gatk4:4.6.1.0--py310hdfd78af_0
SINGULARITY_BIND="-B /lustre10/home,/home"
BCFTOOLS=bcftools

# ワーカー（染色体ごと）のリソース
PARTITION_W=epyc
CPUS_W=8
MEM_W=48G
TIME_W=0-24:00:00
ARRAY_MAXPAR=18   # 同時実行上限（% の数）

# GenomicsDBImport チューニング
BATCH_SIZE=50
READER_THREADS=2

# 結合ジョブのリソース
PARTITION_C=epyc
CPUS_C=2
MEM_C=8G
TIME_C=0-12:00:00

TMPDIR_ROOT=$HOME/working/tmp_gatk
```

ARRAY_MAXPAR … 同時に走らせる染色体ワーカー数の上限。混雑に応じて調整。

CPUS_W と READER_THREADS … ワーカー内の並列のバランス。I/Oが細ければ READER_THREADS を上げ過ぎない。

OUT_ROOT 下から *.g.vcf.gz を探索して SAMPLEMAP（sample名 ↔ gVCFパス） を自動生成します。既存のサンプルマップを使いたい場合は SAMPLEMAP=/path/to/map.txt を渡します。

## 3. ランチャーの動き
### 3.1 事前チェック & 作業ディレクトリ

参照FASTAと染色体リストの存在確認

RUN_DIR に gendb/ vcf/ logs/ を作成

### 3.2 SAMPLEMAP の生成

OUT_ROOT 配下の *.g.vcf.gz を探索

ファイル名から <sample> を切り出し、<sample>\t<絶対パス> の形式で samplemap.txt を作成

>[!TIP]
既にマップがある場合は SAMPLEMAP=/abs/path/to/samplemap.txt を環境変数で渡すとそちらを使用します。

### 3.3 ワーカー（配列ジョブ）の投入
```bash
N=$(grep -Ev '^\s*$|^\s*#' "${CHR_LIST}" | wc -l)
sbatch --array=1-"$N"%${ARRAY_MAXPAR} \
  -J joint_chr -p ${PARTITION_W} --cpus-per-task=${CPUS_W} --mem=${MEM_W} -t ${TIME_W} ...
```

インデックスは 1 始まり（sed -n "${SLURM_ARRAY_TASK_ID}p" で行を取るため）

各タスクは 1 染色体を担当、ログは RUN_DIR/logs/joint_chr_<ジョブID>_<配列ID>.*

### 3.4 全完了後に結合ジョブを依存投入
```bash
sbatch --dependency=afterok:<配列ジョブID> -J joint_concat ...
```
## 4. ワーカーの処理（MODE=worker）
```bash
# 染色体名の取得
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHR_LIST}")"

# 出力先の決定

DBPATH="RUN_DIR/gendb/${CHROM}_db"

VCF_OUT="RUN_DIR/vcf/${CHROM}.vcf.gz"

# 再実行安全
VCF_OUT が既にあれば SKIP して即終了。

# GenomicsDBImport（染色体ごと）

gatk GenomicsDBImport \
  --sample-name-map ${SAMPLEMAP} -L ${CHROM} \
  --genomicsdb-workspace-path ${DBPATH} \
  --batch-size ${BATCH_SIZE} --reader-threads ${READER_THREADS} \
  --tmp-dir ${TMPDIR_ROOT}


# GenotypeGVCFs（染色体ごと）

gatk GenotypeGVCFs \
  -R ${REF_DIR}/${genome}.fa -V gendb://${DBPATH} -L ${CHROM} \
  --tmp-dir ${TMPDIR_ROOT} -O ${VCF_OUT}
tabix -p vcf ${VCF_OUT} || true
```

>[!WARNING]
既存の gendb/${CHROM}_db が残っている状態で再実行すると Import でエラーになる場合があります。
完全な再実行時は該当 gendb/ を消すか、完成済み VCF を残して Import をスキップするロジックを追加する運用も検討ください。

## 5. 結合ジョブの処理（MODE=concat）

染色体ごとの vcf/<chr>.vcf.gz を全て列挙 (vcf_list.txt)
無い染色体があれば エラー終了
```bash
bcftools concat -f vcf_list.txt -O z -o cohort.raw.vcf.gz
tabix -p vcf cohort.raw.vcf.gz
```
>[!TIP]
結合は 連結（concat） であり merge ではありません。同一サンプル集合・非オーバーラップの contig 分割を前提にしています。

## 6. 実行例（投入コマンド）
```bash
# デフォルトで投入
sbatch joint_genotype.sh

# 並列を強め、ワーカーを太くする例
sbatch --export=ALL,ARRAY_MAXPAR=24,CPUS_W=12,MEM_W=64G,READER_THREADS=3 joint_genotype.sh

# 参照・リストを明示
sbatch --export=ALL,REF_DIR=$HOME/ref,CHR_LIST=$HOME/ref/chromosomes.list,OUT_ROOT=$HOME/working/pipeline_out joint_genotype.sh
```
## 7. モニタリング & ログ
```bash
squeue -u $USER | egrep 'joint_chr|joint_concat'
tail -f ${RUN_DIR}/logs/joint_chr_* | sed -u -n 's/.*INFO.*/&/p'
```

各染色体の進捗は RUN_DIR/logs/ に出ます。

## 8. よくある調整ポイント

同時並列：ARRAY_MAXPAR を増減（クラスタ混雑・総メモリ= MEM_W × ARRAY_MAXPAR を意識）

ワーカー内の並列：CPUS_W と READER_THREADS を整合

例）CPUS_W=8 なら READER_THREADS=2〜3 程度から

一時領域：TMPDIR_ROOT を速いローカルストレッチにすると安定

JavaのCPU認識（任意）：
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK}" をセットすると GC スレッドがコア過剰認識しにくくなります。

結合の圧縮並列（任意）：
bcftools concat は -@ で BGZF 圧縮スレッドを増やせます（例：-@ 4）。その場合は CPUS_C も増やす。

## 9. 出力物
```bash
RUN_DIR/
├── gendb/                 # 染色体ごとの GenomicsDB ワークスペース
├── vcf/
│   ├── <chr>.vcf.gz       # 染色体ごとのジェノタイプ結果
│   └── <chr>.vcf.gz.tbi
├── cohort.raw.vcf.gz      # 連結結果
├── cohort.raw.vcf.gz.tbi
└── logs/
    ├── joint_chr_<JID>_<AID>.out/err
    └── joint_concat_<JID>.out/err
```
## 10. 再実行の考え方

染色体 VCF が存在すれば ワーカーは SKIP。不足分だけ再計算されます。

破損した gendb/ があると Import で失敗することがあります。問題の染色体について
RUN_DIR/gendb/<chr>_db を削除して再投下するか、Import スキップのロジック追加を検討してください。
