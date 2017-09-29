% code to test utide

fn = '/Users/pm7/Documents/ptools_data/tide/CO-OPS__9441102__hr.csv';

a = importdata(fn);
t = a.textdata(2:end,1);
dt = datenum(t);
z = a.data(:,1);

coef = ut_solv( dt, z, [], 50, 'auto', 'NoDiagn', 'RunTimeDisp', 'nnn');

[zfit, ~] = ut_reconstr(dt, coef);

dt0 = dt - dt(1);

plot(dt0 , z, '-r', dt0, zfit, '-b', dt0, z-zfit, '-k');

%%
nc = size(coef.name, 1);
for ic = 1:nc
    c = coef.name(ic);
    if strcmp(c,'M2')
        disp(c)
    end        
end


