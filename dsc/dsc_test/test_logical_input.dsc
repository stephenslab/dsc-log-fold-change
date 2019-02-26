# test logical input variable
#
# simple example of input variable can be T or F
# dsc allows this type of input no problem
test: R()
  x: T, F
  $y: 1

DSC:
  run: test
