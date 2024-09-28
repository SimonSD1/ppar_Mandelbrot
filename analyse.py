import matplotlib.pyplot as plt
import numpy as np

cores = np.array([1, 2, 3,4, 5, 6,7,8,9, 10, 20, 30,40, 50, 60, 70, 80, 90, 100])
times = np.array([  62.5103 #1
                      , 34.8349 #2
                      , 29.8332 #3
                      , 23.662  #4
                      , 19.9077 #5
                      , 16.8869 #6
                      , 15.1763 #7
                      , 13.3107 #8
                      , 12.4016 #9
                      , 11.0534 #10
                      , 5.80053 #20
                      , 3.94962 #30
                      , 2.94457 #40
                      , 2.38665  #50
                      , 2.01865  #60
                      , 1.73365  #70
                      , 1.52102  #80
                      , 1.34715  #90
                      , 1.2074   #100
                      ])

plt.figure(figsize=(10, 6))

plt.plot(cores, times, marker='o', linestyle='-', color='b', label="Computing Time")

plt.xlabel('Number of Processors (n)')
plt.ylabel('Total Computing Time (seconds)')
plt.title('MPI Performance: Computing Time vs Number of Processors')

# Optional: Add grid and legend
plt.grid(True)
plt.legend()

# Step 4: Show the plot
plt.savefig("plot.png")


