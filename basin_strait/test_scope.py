"""
Code to test scope.

RESULT: this demonstrates the unexpected behavior that
when a function changes a dict entry it is changed
outside the scope of the function, even if we do not
return the dict!

The reason, apparently, is that:

"the dictionary insertion is not an assigment, but a
method call.  In fact, inserting a key-value pair into
a dictionary is equivalent to calling the __setitem__ method
on the dictionary object."

Interesting!

Even more interesting, or disturbing, is that the function
call h() below also modifies the dict, despite the fact that
I did not pass it anything. h() only does this if we specifically
modify dict "d" and not some other name.  So why does the function
know about the dict object??

Apparently "You don't assign to dictionaryVar, you assign to
dictionaryVar['A'].  So it's never being assigned to, so it's
IMPLICITLY GLOBAL (emphasis mine).  If you were to actyually assign
to dictionaryVar, you'd get the behavior you were 'expecting'."

"""
import numpy as np

def f(a):
    a = 2*a
    return a
    
def g(d):
    for k in d.keys():
        d[k] = d[k]*2
        
def h():
    for k in d.keys():
        d[k] = d[k]+10
    
        
a = 1
b = np.arange(3)
d = {'x':1}

print('before')
print(a)
print(b)
print(d)

f(a)
f(b)
g(d)

print('\nafter without return from f()')
print(a)
print(b)
print(d)

a = f(a)
b = f(b)
g(d)

print('\nafter with return from f()')
print(a)
print(b)
print(d)

h()
print('\nafter call to h()')
print(d)
