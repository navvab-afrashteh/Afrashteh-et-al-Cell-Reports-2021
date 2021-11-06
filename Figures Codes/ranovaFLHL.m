function prT = ranovaFLHL(annovaVar)
% run repeated measure anova on data from FL or HL stim

animals = {'animal 1','animal 2','animal 3','animal 4','animal 5'}';
subnetworks = {'SCC', 'L', 'M', 'H'};
t = table(animals,annovaVar(:,1),annovaVar(:,2),annovaVar(:,3),annovaVar(:,4),...
'VariableNames',{'animals','meas1','meas2','meas3','meas4'});
Meas = table(subnetworks','VariableNames',{'subnetworks'});
rm = fitrm(t,'meas1-meas4~1','WithinDesign',Meas);
ranovatbl = ranova(rm);
% estimated partial effec size (etta^2)
etta2 = ranovatbl.SumSq(1)/(ranovatbl.SumSq(2)+ranovatbl.SumSq(1));
% check for sphericity. if pvalue<0.05 then use correction. For correction
% calculate epsilon, depending on this value use proper correction type.
sph = mauchly(rm);
fprintf('ranova: F = %f\n',ranovatbl.F(1))
fprintf('ranova: etta2 = %f\n',etta2)
if sph.pValue>=0.05
    pValue = ranovatbl.pValue(1);
else
    ep = epsilon(rm);
    if ep.GreenhouseGeisser>0.75
        pValue = ranovatbl.pValueHF(1);
    else
        pValue = ranovatbl.pValueGG(1);
    end
end
fprintf('ranova: pValue = %0.10f\n',pValue)
post_hoc = 'hsd';
tbl = multcompare(rm,'subnetworks','ComparisonType', post_hoc);
prT = tbl.pValue([11,12,10,5,4,7]);
if (pValue>=0.05)
    prT(:) = 1;
    tbl.pValue([11,12,10,5,4,7]) = 1;
end

tbl([11,12,10,5,4,7],{'subnetworks_1','subnetworks_2','pValue'})
