# test logical input variable
#
# simple example of input variable can be T or F
# dsc allows this type of input no problem
data_poisthin: R()
  seed: R{2:3}
  nsamp: 90
  ngene: 1000
  prop_null: .5, .9, 1
  shuffle_sample: T, F
  gselect: "random"
  signal_dist: "bignormal"
  $Y1: 1
  $Y2: 1


DSC:
  define:
    data: data_poisthin
  run:
    data
