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

When succefully executed, you should see something similar to the following:
```sh
Distance between  NZ_AKBW01000001.1  and  .14412_3_82.9 is:  0.9994697773064687
Distance between  NZ_AKBW01000001.1  and  .14412_3_82.10 is:  0.9994666666666666
Distance between  NZ_AKBW01000001.1  and  .14412_3_82.11 is:  1.0
Distance between  NZ_AKBW01000001.1  and  .14412_3_82.12 is:  0.9994706193753309
Distance between  NZ_AKBW01000001.1  and  .14412_3_82.13 is:  1.0
Distance between  NZ_AKBW01000001.1  and  AE007317.1 is:  0.26956521739130435
Distance between  .14412_3_84.1  and  .14412_3_84.2 is:  0.993798449612403
Distance between  .14412_3_84.1  and  .14412_3_84.3 is:  0.9943298969072165
Distance between  .14412_3_84.1  and  .14412_3_84.4 is:  0.994335736354274
Distance between  .14412_3_84.1  and  .14412_3_84.5 is:  0.9948293691830403
Distance between  .14412_3_84.1  and  .14412_3_84.6 is:  0.9974437627811861
Distance between  .14412_3_84.1  and  .14412_3_84.7 is:  0.9974385245901639
Distance between  .14412_3_84.1  and  .14412_3_84.8 is:  0.9969056214543579
Distance between  .14412_3_84.1  and  .14412_3_84.9 is:  0.9989801121876594
...
```

## Effect of Sketch Size
### Effect on distance
We tried to analyse the effect of the sketch size on the calcualted distance. So, we took two random sequences `.14412_3_84.17` and `.14412_3_82.13` and calulcated their distances using different values of the sketch size, and we got the following results:

![sketch_size_to_distance](https://github.com/user-attachments/assets/a23eecc0-a3db-47ac-ad46-fc9a74aa6324)

We can see that for sketch sizes less than **1000**, the distance was unpredictable. But starting from **1000**, the distance remains the same at `0.08561643835616439`.
### Effect on calulation time
The following line plot illustrates the effect of the sketch size on execution time:

![execution_time](https://github.com/user-attachments/assets/b0905776-85eb-4daf-aae0-fda038c3caa0)
