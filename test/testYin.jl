using PhiVoc.Yin

x = zeros(16)
x[[1,9]] .= 1

max_lag = 4

# difference function
@test length(PhiVoc.Yin.AMDF(x, max_lag)) ==  max_lag 
# minimum finder
y = [1.,.5,.2,.1,.05,.01,.05]
@test PhiVoc.Yin.findFirstMinimum(y, .1) == 6
# paraboilc interpolation
@test PhiVoc.Yin.parabolicInterp([1,0,1],2) == (2.0, 0.0)
@test PhiVoc.Yin.parabolicInterp([5,1,1],2) == (2.5, 0.5)
@test PhiVoc.Yin.parabolicInterp([0.,2,4],1) == (1.0, 0.0)

