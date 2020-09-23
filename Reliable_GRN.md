---
title: "# 酵母の信頼できる遺伝子ネットワーク"
tags: ""
---

Yeastract <http://www.yeastract.com/> のRegulation matrixからダウンロードした。

-   TF-consensus list にある全ての転写転写因子の名前をRegulation matrix のTranscription factorsの欄にペーストする。

- **both.csv:  活性化と抑制の両方
- **activator.csv: 活性のみ
- **inhibitor.csv: 抑制のみ
パラメータ

**both.csv  

> -   Documented
> -   TF acting as activator
> -   TF acting as represser
> -   DNA binding plus expression evidence
> -   Unstressed group condition
> -   Consider User inserted TF list

**activator.csv  

> -   Documented
> -   TF acting as activator
> -   DNA binding plus expression evidence
> -   Unstressed group condition
> -   Consider User inserted TF list

**inhibitor.csv  

> -   Documented
> -   TF acting as represser
> -   DNA binding plus expression evidence
> -   Unstressed group condition
> -   Consider User inserted TF list



より厳しい(確実な)条件でフィルターしたデータ
RegulationTwoColumnTable_both_strict.tsv

> -   Documented
> -   TF acting as represser
> -   DNA binding and expression evidence
> -   Unstressed group condition
> -   Consider User inserted TF list 


Gene -> ORF の変換に関して

下記コマンドでカンマを除きコピーした後
cut -d ',' -f 2 ./YRGRNlist.csv  | sed -e 's/"//g' | pbcopy

yeastract に突っ込む
