from random import random
with open ("dots.dat", "w") as file:
    for i in range (1000) :
        print (i , (3 * i**3 - 2 * i**2 + 7 * i - 9) +( random () - 0.5) *10 , file= file )