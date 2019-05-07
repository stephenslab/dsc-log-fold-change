# inherent arguments from previous modules
#
test: R()
  x: T, F
  x1: 1,2
  x2: 3,4
  $y: 1

#test2(test):
#  $y: 2

test2(test):

DSC:
  run: test2
  output: test_def_modules

#exec_path: /project2/gilad/joycehsiao/dsc-log-fold-change/dsc/modules
