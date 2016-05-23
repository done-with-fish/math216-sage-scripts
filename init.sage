import getpass

latex.matrix_delimiters(left='[', right=']')

def my_test():
    print "hello math 216!"
    return None

if getpass.getuser() == 'brian':
    load('/home/brian/ownCloud/code/math216-sage-scripts/rref-steps.sage')
else:
    load('/home/math216/.scripts/math216-sage-scripts/rref-steps.sage')
