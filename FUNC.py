#################################################
## 3 # HELPER FUNCTIONS, CLASSES AND SHORTCUTS ##  -> @FUNC <-
#################################################

import math

#----+------------------+
## A | STRING FUNCTIONS |
#----+------------------+


# Split a string
def spl(x):
    return x.split()


# Split each argument in a list
def nsplit(*x):
    return [i.split() for i in x]


# Make a dictionary from two lists
def hash(x, y):
    return dict(list(zip(x, y)))


# Function to reformat pattern strings
def pat(x, c="."):
    return x.replace(c, "\x00").split()


# Function to generate formatted strings according to the argument type
def formatString(i):
    if type(i) == str:
        return i
    if type(i) == int:
        return "%5d" % i
    if type(i) == float:
        return "%8.5f" % i
    else:
        return str(i)


#----+----------------+
## B | MATH FUNCTIONS |
#----+----------------+


def cos_angle(a, b):
    p = sum([i*j for i, j in zip(a, b)])
    q = math.sqrt(sum([i*i for i in a])*sum([j*j for j in b]))
    return min(max(-1, p/q), 1)


def norm2(a):
    return sum([i*i for i in a])


def norm(a):
    return math.sqrt(norm2(a))


def distance2(a, b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
