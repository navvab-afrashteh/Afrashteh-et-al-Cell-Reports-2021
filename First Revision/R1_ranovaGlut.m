function prT = R1_ranovaGlut(annovaVar)
% run repeated measure anova on data from FL or HL stim

animals = {'animal 1','animal 2','animal 3','animal 4','animal 5'}';
subnetworks = {'S', 'SCC', 'L', 'H'};
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
    fprintf('ranova: pValue = %0.10f\n',ranovatbl.pValue(1))
else
    ep = epsilon(rm);
    if ep.GreenhouseGeisser>0.75
        fprintf('ranova: pValueHF = %0.10f\n',ranovatbl.pValueHF(1))
    else
        fprintf('ranova: pValueGG = %0.10f\n',ranovatbl.pValueGG(1))
    end
end
post_hoc = 'hsd';
tbl = multcompare(rm,'subnetworks','ComparisonType', post_hoc)
%prT = tbl.pValue([6,5,3]);
prT = tbl.pValue([9,8,7,11,10,1]);

tbl([9,8,7,11,10,1],{'subnetworks_1','subnetworks_2','pValue'})
