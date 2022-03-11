clear

I_s = 0.01 * 10 ^ -12; % A
I_b = 0.1 * 10 ^ -12; % A
V_b = 1.3; % V
G_p = 0.1; % 1/Ohm

noise = 0.20; % PErcent varriation of current

V = linspace(-1.95, 0.7, 200);

% Ideal Current 
I = I_s*(exp(1.2/0.025*V)-1) - I_b*(exp(-1.2/0.025*(V+V_b))-1)+ G_p * V ;

% Current with noise
I_rand = zeros(size(I));
randNum = -noise + 2*noise * rand(200, 1); %Random numbers between -0.2 and +0.2
for n = 1:200
 I_rand(n) = I(n) + I(n) * randNum(n);
end


% Curve Fitting
x = linspace(-1.95,0.7);
fit_I_rand_4 = polyfit(V,I_rand,4);
fit_I_rand_8 = polyfit(V,I_rand,8);
fit_I_4 = polyfit(V,I,4);
fit_I_8 = polyfit(V,I,8);

f0 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - 0.1 * 10 ^ -12*(exp(1.2*(-(x+1.3))/25e-3)-1)');
f1 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');


ff0 = fit(V',I',f0, 'StartPoint', [ I_s/2, G_p*3/2]);
If0 = ff0(x);
ff1 = fit(V',I',f1, 'StartPoint', [I_s/2 G_p*3/2 I_b/2 V_b*3/2]);
If1 = ff1(x);

% Plotting 
tiledlayout(2,1)
nexttile
plot(V, I, 'marker', '.');
hold on 
plot(x, If0)
title("Diode Current, I_s and G_p Parameters Optimized" )
xlabel('V_d (V)');
ylabel('I_d (A)');
nexttile
plot(V, I, 'marker', '.');
hold on 
plot(x, If1)
title("Diode Current, all Parameters Optimized" )
xlabel('V_d (V)');
ylabel('I_d (A)');
hold off

figure
tiledlayout(2,2)
nexttile
plot(V, I);
hold on
plot(x, polyval(fit_I_4, x))
hold on
plot(x, polyval(fit_I_8, x))
xlabel('V_d (V)');
ylabel('I_d (A)');

nexttile
semilogy(V, I);
xlabel('V_d (V)');
ylabel('log(I_d) (log(A))');

nexttile
plot(V, I_rand);
hold on
plot(x, polyval(fit_I_rand_4, x))
hold on
plot(x, polyval(fit_I_rand_8, x))
xlabel('V_d (V)');
ylabel('I_d (A)');

nexttile
semilogy(V, I_rand);
xlabel('V_d (V)');
ylabel('log(I_d) (log(A))');




inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs)
view(net)
Inn = outputs