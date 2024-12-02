# genome-distance
Efficient way of calculating the distance between two genome sequences using the Jaccard Distance 

## Run
To the run the code, follow the below steps:
```sh
git clone https://github.com/waterflow80/genome-distance.git
cd genome-distance
python3.11 -m venv venv
./venv/bin/activate
pip install -r requirements
python3.11 genome_distance.py
# You might need to include the python interpreter of the venv
# ./venv/bin/python3 genome_distance.py
```

When succefully executed, you should see the following output:
```sh
Distance between  TIGR4.fa  and  14412_3#84.contigs_velvet.fa is:  0.9884805185710432
Distance between  TIGR4.fa  and  14412_3#82.contigs_velvet.fa is:  0.9884531595043302
Distance between  TIGR4.fa  and  R6.fa is:  0.9743687595195952
Distance between  14412_3#84.contigs_velvet.fa  and  14412_3#82.contigs_velvet.fa is:  0.5608253422006775
Distance between  14412_3#84.contigs_velvet.fa  and  R6.fa is:  0.9887006693223716
Distance between  14412_3#82.contigs_velvet.fa  and  R6.fa is:  0.9886752171486473
```
