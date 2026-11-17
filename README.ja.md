# ToMMo NMR Quantification

NMRのBINスペクトルから**血漿中50代謝物**の濃度を推定し、µM単位の濃度表を作成するパイプラインです。

**BINスペクトル → 化合物別の特徴量 → 濃度推定 → 指定バッチの補正 → 10SD除外 → 低値除外**

[実行手順](docs/setup.md) · [全50項目のモデル一覧](config/model_columns.tsv) · [特徴量・バッチ補正](docs/curation.md)

## 定量方法

入力は、検体IDと共通の[ppm軸](reference/BinTable_Header.tsv)に対応する37,000列のスペクトル強度です。モデルと同じppm軸・強度スケールのスペクトルを使用します。

化合物ごとに指定した領域を抽出し、対応するモデルで濃度を推定します。

- **46化合物は線形回帰**：BIN強度の総和、ピーク高さ、または補正したピーク高さを特徴量として使用します。
- **4化合物は多層パーセプトロン（MLP）**：メチオニン、グリセロール、バリン、アラニンについて、指定領域のスペクトル強度を入力します。

線形回帰では `濃度 = 係数 × 特徴量 + 切片` として計算します。BIN総和にはppm幅を追加で乗じません。必要な化合物では、参照シグナル、ピーク位置の変動、重なりに対応した特徴量処理を行います。推定時は同梱の係数・重みを適用します。

[モデル一覧](config/model_columns.tsv)に、全50化合物の方法、入力列、ppm、特徴量、係数または重みファイルを記載しています。TSV列番号は**IDを第1列として数え、領域の両端を含みます**。実行用の定義は [models.json](config/models.json) に格納しています。

## 濃度表の整備

指定したバッチには、領域変更やメディアン差に基づく補正を適用できます。検体とバッチの対応・設定は入力ファイルで指定します。

必要な濃度補正の後、次の順序で処理します。

1. 化合物別の基準分布の**平均±10SDの外側**をNAにします。境界と等しい値は残します。
2. 残った値のうち、化合物別の[運用上の下限値以下](config/reporting_cutoffs.json)をNAにします。

平均・SDには同梱の[基準統計表](config/reference_stats.tsv)を用い、入力検体から計算し直しません。各段階の濃度表と除外理由を保存します。

## 実行

```bash
conda env create -f environment.yml
conda activate tommo-nmr-quantification
bash scripts/run_synthetic_demo.sh
```

人工スペクトルによる実行例です。実データの入力形式とバッチ設定は[実行手順](docs/setup.md)を参照してください。計算にはfloat64を用い、閲覧用の表だけを丸めます。

引用情報は [CITATION.cff](CITATION.cff) に記載しています。

公開コードと同梱モデルの係数・重みには [MIT License](LICENSE) を適用します。著作権表示は ToMMo NMR Quantification Project Contributors です。コホート参加者データの利用条件は別途適用されます。
