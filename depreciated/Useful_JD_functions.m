%% Useful functions to review

getrow
pivottable
tapply
lineplot
scatterplot
boxplot
traceplot

e.g.
pivottable(D.BN,D.digit,D.EVal,'mean');
pivottable(D.BN,D.digit,D.digit,'length(x)');
D
T = tapply(D,{'BN'},{'EVal','mean'})
T = tapply(D,{'BN','digit'},{'EVal','mean'})

lineplot(T.BN,T.EVal)
lineplot(T.BN,T.EVal,'split',T.digit)
lineplot(T.BN,T.EVal,'split',T.digit,'style_thickline')
lineplot(T.BN,T.EVal,'split',T.digit,'style_thickline','leg','auto')
help lineplot