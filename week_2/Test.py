import math

# veerconstante, k (N/m)
k = 40

# massa, m (kg)
m = 1.845 * 10 ** -6

# gamma, c (kg/s)
c = 0.042953

b = 2 * c * m

result = math.sqrt(4*m*k)

if b < result:
    print("Ondergedempt")
elif b == result :
    print("Kritisch gedempt")
elif b > result:
    print("Overgedempt")
    
print("result")