# Format files for IPA
echo $'ID\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB' > D21_WT.txt
tail -n +2 ../01.Differential_Expression/D21_WT.xls | sed 's/|/\t/g' | cut -d $'\t' -f2,4,5,6,7,8,9 >> D21_WT.txt
echo $'ID\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB' > D21_D14.txt
tail -n +2 ../01.Differential_Expression/D21_D14.xls | sed 's/|/\t/g' | cut -d $'\t' -f2,4,5,6,7,8,9 >> D21_D14.txt
echo $'ID\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB' > D14_WT.txt
tail -n +2 ../01.Differential_Expression/D14_WT.xls | sed 's/|/\t/g' | cut -d $'\t' -f2,4,5,6,7,8,9 >> D14_WT.txt
