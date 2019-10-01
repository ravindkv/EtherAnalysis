inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
cd 13TeV
mkdir $outAnalDir
cd ../

echo "INPUT FILE: " $inNtupleFile
cp Analyzer.C ../
sed -i "s:inputFile:$inNtupleFile:g" Analyzer.C
sed -i "s:outputFile:$outAnalFile:g" Analyzer.C
sed -i "s:outputDir:$outAnalDir:g" Analyzer.C

root -l -b <<EOF
.L Analyzer.C+
Analyzer t
t.processEvents()
.q
EOF
rm *.so
rm *.d
mv ../Analyzer.C .
