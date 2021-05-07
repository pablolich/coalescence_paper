monod = function(a, N, k){
  return(a*N/(k + N))
}
biochemistry = function(x, b, n){
  return(x^n/(b^n + x^n))
}
x = seq(100)
b = 2
n = 1
bio1 = biochemistry(x, b, n)
n = 2
bio2 = biochemistry(x, b, n)
n = 3
bio3 = biochemistry(x, b, n)
plot(x, bio1, type = 'l')
lines(x, bio2, col = 'red')
lines(x, bio3, col = 'green')

growth = function(t, N0, r){
  return(N0*exp(r*t))
}

N0 = 100
t = seq()
r = 2
growth1 = growth(t, N0, r)
r = 3
growth2 = growth(t, N0, r)
plot(t, growth2, type = 'l')
lines(t, growth1, col = 'green')

fish = function(Linf, k, x){
  return(Linf*(1-exp(-k*x)))
}

Linf = 20
k = 1
x = seq(100)
fish1 = fish(Linf, k, x)
k = 0.1
fish2 = fish(Linf, k, x)
plot(x, fish1, type = 'l')
lines(x, fish2, col = 'red')

par(mfrow = c(3,1))
x = linspace(0.1,9, 100)
plot(x, exp(x)/(x^2), type = 'l', lwd = 5, ylab = 'Ratio', cex.lab = 3, cex.axis = 2)
legend(x = 2, y = 80, legend=c(expression(paste(e^x/x^2))),  pch = 20, cex=2)
x = linspace(0.1,65, 150)
plot(x, exp(x)/(x^10), type = 'l', lwd = 5, ylab = 'Ratio', cex.lab = 3, cex.axis = 2)
legend(x = 10, y = 1e10, legend=c(expression(paste(e^x/x^10))), pch = 20, cex=2)
x = linspace(0.1,10, 1500)
plot(x, log10(x)/x, type = 'l', lwd = 5, ylab = 'Ratio', cex.lab = 3, cex.axis = 2)
legend(x = 2, y = -4, legend=c(expression(paste(log(x)/x))), pch = 20, cex=2)



fac = function(l){
  return( l/(1-l*(1-l)) - l )
}

comp = function(l){
  return( (1-l)/(1-l*(1-l)) )
}

new_fac = function(l){
  return( -1/(log(l*(1-l))) )
}

new_comp = function(l){
  return((1-l) - 1/(log(l*(1-l))) )
}

new_fac2 = function(l){
  return( -l/(log(l*(1-l))) )
}

new_fac = function(l){
  return((l-1)*l/(log(l*(1-l))) )
}

new_comp = function(l){
  return( 1 - l + ((l-1)*l/(log(l*(1-l)))) )
}

new_fun = function(l){
  return((l-1)/log(l))
}



plot(l, new_comp(l), type = 'l', lty = 'dashed', col = 'black', lwd = 4)
lines(l, new_fac(l), type = 'l', col = 'green',lwd = 2)
lines(l, 1-new_fun(l), col = 'grey', lwd = 2)
  3 vqqqlines(l, new_fun(l), type = 'l', col = 'red', lwd = 2)
abline(v = 2 - sqrt(3), col = 'black', lty = 'dashed')
abline(v = 0.35, lty = 'dashed')
plot(l, fac(l)- comp(l), type = 'l')


l = seq(from=0, to = 0.987, by = 0.001)
df = data.frame(l = l, B = new_fun(l), Bi = 1-new_fun(l))
mdf = melt(df, measure.vars = c('B', 'Bi'))
ggplot(data = mdf)+
  geom_line(aes(x = l, y = value, 
                group = variable,
                color = variable),
            size = 1.5)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.5),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 20),
        legend.text.align = 0,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 10),
        aspect.ratio = 1)+
  scale_x_continuous(expand = c(0.007, 0.007)) +
  scale_y_continuous(expand = c(0.007, 0.007)) +
  scale_color_manual(values = c('black', 'grey'), 
                     labels = c(expression(B), expression(1-B)))+ 
  geom_point(aes(x = 1, y = 1),
             shape = 1, size = 2.5)+
  geom_point(aes(x = 1, y = 0),
             shape = 16, size = 2.5)+
  geom_segment(x = 1, xend = 1, y = 0, yend = 1, linetype = 'dashed')+
  geom_segment(x = 0, xend = 1, y = 1, yend = 0, linetype = 'dashed')+
  labs(x = expression(l), y = expression(f(l)))+
  ggsave(filename = '../sandbox/biotic_factor.pdf',
         height = 5, 
         width = 5)
