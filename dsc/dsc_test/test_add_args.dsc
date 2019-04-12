# see if dsc query would change after adding
# additional input variables
#
# simple example of input variable can be T or F
# dsc allows this type of input no problem
test: R()
  x: T, F
  x1: 1,2
  x2: 3,4
  $y: 1

DSC:
  run: test
  output: test_add_args

