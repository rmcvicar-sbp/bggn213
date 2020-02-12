 data.frame(id = letters[1:10], x = 1:10, y = 11:20)
data.frame()
data.frame(x=1:10, y=11:20)
data.frame(letters[1:10], x=1:10, y=11:20)
data.frame(id = letters[1:10], x = 1:10, y = 11:20)



# Code the control-flow construct
if (li >= 15 & fb >= 15) {
  sms <- 2 * (li + fb)
} else if (li < 10 & fb < 10) {
  sms <- 0.5 * (li + fb)
} else {
  sms <- li + fb
}
