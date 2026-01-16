import numpy
import matplotlib.pyplot as plt

data = {6:60, 7:75, 8:95, 9:116, 10: 140, 11:165, 12:196, 13:225, 14:270, 15:290, 16:335, 17:380, 18:420, 19:465, 20:515}
x = []
y = []
for k in data.keys():
    x.append(k)
    y.append(data[k])

mymodel = numpy.poly1d(numpy.polyfit(x, y, 2))
print(mymodel)

myline = numpy.linspace(1, 22, 100)

plt.scatter(x, y)
plt.plot(myline, mymodel(myline))
plt.show()
