#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

#I think I can do this w/ a for loop. I'll need it to do 4 things:
#establish a range between 0 and 99
#make if-statement to check if remainder obtained when the number divided by 7 is one or not.
#If the remainder isn’t equal to 0, the number isn’t divisible by 7.
#only numbers satisfying the condition are printed

for i in range(0,99): 
   if(i%7!=0): 
       print(i)
