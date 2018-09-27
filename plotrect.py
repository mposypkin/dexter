import matplotlib.pyplot as plt
import sys
import re

plt.axes()

# pat = re.compile('\[(\-|\d|\.)*\:(\-|\d|\.)*\]')
pat = re.compile('\[(?:\d|-|\.)*\:(?:\d|-|\.)*\]')
for line in sys.stdin:
    # print(line)
    if re.match("^#.*", line):
        continue
    arr =pat.findall(line)
    # print(arr)
    # ends = re.split('\[|\:|\]', I)
    # I.strip('[')
    I = arr[0]
    l = len(I)
    ends = re.split('\:', I[1:l-1])
    # print(ends)
    a1 = float(ends[0])
    b1 = float(ends[1])

    I = arr[1]
    l = len(I)
    ends = re.split('\:', I[1:l-1])
    # print(ends)
    a2 = float(ends[0])
    b2 = float(ends[1])


    rectangle = plt.Rectangle((a1, a2), b1 - a1, b2 - a2, fc='y', fill = True, color='red')
    plt.gca().add_patch(rectangle)


#circle = plt.Circle((0, 0), radius=3, fc='y', fill = False)
#plt.gca().add_patch(circle)
plt.axis('scaled')
plt.show()
