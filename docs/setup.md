# 実行手順

リポジトリのルート（README.mdがある場所）で実行します。

## 環境と動作確認

```bash
conda env create -f environment.yml
conda activate tommo-nmr-quantification
python scripts/verify_release.py
python -m unittest discover -s tests -v
bash scripts/run_synthetic_demo.sh
```

指定環境はPython 3.9.12、NumPy 1.21.5、pandas 1.4.2、SciPy 1.9.3です。Python環境を既に用意している場合は `python -m pip install -r requirements.txt` でも導入できます。人工6検体の例が全50項目と両方の除外段階を実行し、最後に `PASS ... samples=6 metabolites=50` を表示します。

オフライン導入では、同じOS・CPU・Pythonの接続可能な環境でwheelを用意して持ち込みます。

```bash
python -m pip download --only-binary=:all: --dest wheelhouse -r requirements.txt
# 持ち込み先のPython環境で実行
python -m pip install --no-index --find-links wheelhouse -r requirements.txt
```

## 入力

* **BINファイル**：タブ区切り、ID＋37,000列。共通の `reference/BinTable_Header.tsv` と同じ軸・強度スケールを使います。`--format header` は1行目のppmを照合、`raw8` は先頭8行を読み飛ばし、`none` はヘッダーなしとして読みます。
* **基準分布**：同梱の `config/reference_stats.tsv` を既定で使用します。全50項目の平均・標本SDをµM、小数2桁で保持しています。別の統計表を指定する場合のみ、下の3列のTSVを用意して `--reference-stats` に渡します。

```text
TCN	mean_uM	sd_uM
```

* **バッチ対応表**：`SampleID` と `BatchID` の2列のTSV。
* **バッチ設定**：各バッチのN,N-Dimethylglycine（TCN000042）の基準ピーク位置。人工データ生成後の `runs/synthetic_demo/input/groups.json` に設定形式の例があります。

`tcn42_reference_column` は、入力TSVの列を0から数えた位置（**ID列=0、最初のBIN=1**）です。モデル一覧の1始まりのTSV列番号より1小さい値を設定します。基準位置は各バッチのスペクトルを目視確認して決めます。以降のピークパターンの判定、特徴量調整、濃度計算は自動です。

モデル一覧のBIN領域は、コードではIDを除くスペクトルの0始まりインデックスで保持します。対応式は `TSV列番号 = BINインデックス + 2` です。

## 一括実行

```bash
python -m nmr_cohort run \
  --bin data/spectra.tsv --format header \
  --batch-map data/batch_map.tsv --groups data/groups.json \
  --out runs/cohort

python scripts/check_run.py runs/cohort/quantified runs/cohort/filtered
```

1バッチだけなら、`--batch-map` と `--groups` の代わりに `--reference-column` で基準位置を直接指定できます。出力先は既存フォルダを上書きしません。

| 出力 | 内容 |
|---|---|
| `quantified/raw.tsv` | 固定モデルの濃度推定値、50列 |
| `quantified/features.tsv` | スカラー特徴量、46列 |
| `quantified/tcn42_diagnostics.tsv` | TCN000042の判定モード・処理状態 |
| `filtered/after_10sd.tsv` | 基準分布の10SD除外後 |
| `filtered/concentrations.tsv` | 続いて低値除外を適用した最終表 |
| `filtered/display_2dp.tsv` | 最終表の閲覧用、小数2桁 |
| `filtered/mask_flags.tsv` | 除外段階と理由 |
| 各段階の `run.json` | 完了状態、件数、入力・出力のハッシュ |

入力IDは既定ではそのまま保持します。ピークフィッティングの数値エラーをNAとして継続する場合だけ `--on-fit-error na` を追加します。

## 指定バッチへの補正

固定済みのメディアン補正量を適用するには、一括実行に `--batch-parameters data/batch_parameters.json` を追加します。補正は10SD除外より前に行います。補正量の作成コマンドは次のとおりです。

```bash
python -m nmr_cohort calibrate-batches \
  --concentrations data/calibration_concentrations.tsv \
  --batch-map data/calibration_batch_map.tsv \
  --rules data/batch_rules.json --out data/batch_parameters.json
```

ルール形式は `examples/batch_rules.json` を参照してください。領域変更はバッチ設定の `roi_overrides` に記述できます（例：`"TCN000028": {"roi_ppm": {"high": 3.957, "low": 3.948}, "selector": "bracketed"}`）。通常は全検体に共通のモデル定義を用い、設定したグループだけを変更します。

## ファイルを編集した場合

コード・README・引用情報などを更新した後は、`python scripts/update_manifest.py` でチェックサムを更新し、変更したファイルと `MANIFEST.json` を一緒にコミットしてください。続いて `python scripts/verify_release.py` で確認できます。
