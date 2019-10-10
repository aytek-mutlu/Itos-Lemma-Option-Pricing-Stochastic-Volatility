%%%Ito's lemma

S0=1;
r=0;
sigma = 0.15;
t0=0;
t_N = 1;
NumSteps = 252;
NumPaths = 5000;
S = zeros(NumPaths,NumSteps);
s = zeros(NumPaths,NumSteps);
S(:,1)=1;   
s(:,1)=log(1);
time_step = (t_N-t0)/NumSteps;

for i=2:NumSteps
    brownian = randn(NumPaths,1);
    S(:,i) = S(:,i-1).*(1 + r*time_step + sigma*brownian*sqrt(time_step));
    s(:,i) = s(:,i-1) + (r-power(sigma,2)/2)*time_step + sigma*brownian*sqrt(time_step);
end

S_hat = exp(s);

S_N = S(:,NumSteps);
S_hat_N = S_hat(:,NumSteps);

S_0 = S(:,1);
S_hat_0 = S_hat(:,1);

%1
subplot(1,2,1)
hist(S_N,500);
title('Distribution of S_N')
subplot(1,2,2)
hist(S_hat_N,500);
title('Distribution of Shat_N')

mean(S_N)
mean(S_hat_N)

std(S_N)
std(S_hat_N)

skewness(S_N)
skewness(S_hat_N)

kurtosis(S_N)
kurtosis(S_hat_N)

%test if same distribution
kruskalwallis([S_N,S_hat_N])

%2
log_S_N_0 = log(S_N) - log(S_0);
log_S_hat_N_0 = log(S_hat_N) - log(S_hat_0);

subplot(1,2,1)
hist(log_S_N_0,500);
title('Distribution of log return of S_N')
subplot(1,2,2)
hist(log_S_hat_N_0,500);
title('Distribution of log return of Shat_N')

mean(log_S_N_0)
mean(log_S_hat_N_0)
%mean should be (r-1/2 sigma^2)
r - 0.5*sigma*sigma

std(log_S_N_0)
std(log_S_hat_N_0)
%std should be sigma
sigma

skewness(log_S_N_0)
skewness(log_S_hat_N_0)

kurtosis(log_S_N_0)
kurtosis(log_S_hat_N_0)

%3
%log returns
logS = log(S);
logShat = log(S_hat);

logreturnS = zeros(NumPaths,NumSteps);
logreturnShat = zeros(NumPaths,NumSteps);



for i=1:NumSteps
    logreturnS(:,i) = logS(:,i) - log(S_0);
    logreturnShat(:,i) = logShat(:,i) - log(S_hat_0);
end


plot(logreturnS(1,:)')
hold on

plot(logreturnShat(1,:)')
legend({'log return of S_N','log return of Shat_N'})

plot(logreturnShat(1,:)' - logreturnS(1,:)')
title('Difference')

%%%%monte carlo
int_rate=0.01;
last_price = 2978.4;
compound_freq = 0.25;
option_maturity = 0.25;
annual_simple_int_rate = power((1+int_rate*compound_freq),1/compound_freq)-1;
sigma = 0.1424;
strike = 3000;
M = 100000;

discount  = exp(-annual_simple_int_rate*option_maturity);
MC = last_price*exp((annual_simple_int_rate - 0.5*power(sigma,2))*option_maturity+sigma*sqrt(option_maturity)*randn(M,1));
payoffs =  discount * max(MC-strike,0);
Ms= [100,500,1000,5000,10000,30000,50000,70000,100000];
payoff_list = [mean(payoffs(1:100,:)),mean(payoffs(1:500,:)),mean(payoffs(1:1000,:)),mean(payoffs(1:5000,:)),mean(payoffs(1:10000,:)),mean(payoffs(1:30000,:)),mean(payoffs(1:50000,:)),mean(payoffs(1:70000,:)),mean(payoffs(1:100000,:))];
[call_price,call_delta] = bs_call(last_price,strike,annual_simple_int_rate,option_maturity,sigma,1);

plot(Ms,payoff_list,Ms,call_price*ones(1,length(Ms)))
legend({'Monte Carlo Simulation','Black-Scholes Price'})


%Dynamic Hedge
strike  = 3000;
spot = 2950;
T = 0.25;
%sigma = 0.1424;
sigma = 0.1924;

int_rate= 0.01;
size = -100*10;
NumScenarios = 5000;
%NumPeriods = 13;
NumPeriods = 252/4;

[call_price,call_delta] = bs_call(spot,strike,int_rate,T,sigma,size);


SpotPrices  = zeros(NumScenarios,NumPeriods+1);
SpotPrices(:,1) = spot;

BSCallPrices = zeros(NumScenarios,NumPeriods+1);
BSCallPrices(:,1) = call_price;

Deltas = zeros(NumScenarios,NumPeriods+1);
Deltas(:,1) = call_delta;

PnL = zeros(NumScenarios,NumPeriods+1);

for i=2:(NumPeriods+1)
    time_to_mat = T*((NumPeriods+1-i)/NumPeriods);
    passed_time = T - time_to_mat;
    SpotPrices(:,i) = SpotPrices(:,i-1) .* exp((int_rate - 0.5*power(sigma,2))*T/NumPeriods+sigma*sqrt(T/NumPeriods)*randn(NumScenarios,1));
    [BSCallPrices(:,i),Deltas(:,i)] = bs_call(SpotPrices(:,i),strike,int_rate,time_to_mat,sigma,size);
    PnL(:,i) = PnL(:,i-1) + exp(int_rate*time_to_mat) * Deltas(:,i-1) .* (SpotPrices(:,i) - SpotPrices(:,i-1));
end

Final_PnL = -exp(int_rate*T) * BSCallPrices(:,1) + PnL(:,NumPeriods+1) + BSCallPrices(:,NumPeriods+1);

plot(BSCallPrices(1,:))
hold on
plot(PnL(1,:))
legend({'Short Call P&L','Delta Hedge P&L'})

hist(Final_PnL,100)
mean(Final_PnL)
std(Final_PnL)

%heston
int_rate=0.01;
last_price = 2978.4;
compound_freq = 0.25;
option_maturity = 0.25;
rate = power((1+int_rate*compound_freq),1/compound_freq)-1;
M = 100000;
delta_t = 1/252;
num_periods = option_maturity / delta_t;


sigma_initial = 0.1424;
strike = 3000;
k = 3;
vol_vol = 0.2;
p=0.7;
long_term_sigma = 0.15;
long_term_var = power(long_term_sigma,2);

stock_prices = zeros(M,num_periods+1);
stock_prices(:,1) = last_price;

variances = zeros(M,num_periods+1);
variances(:,1) = power(sigma_initial,2);

for i=2:(num_periods+1)
    
    g_s = randn(M,1);
    g_v = randn(M,1);
    rand1 = g_s;
    rand2 = p*g_s + sqrt(1-power(p,2))*g_v;
    
    stock_prices(:,i) = stock_prices(:,i-1).*(1 + rate*delta_t + sqrt(variances(:,i-1)*delta_t).*rand1);
    variances(:,i) = variances(:,i-1) + k*(long_term_var - variances(:,i-1))*delta_t + vol_vol*sqrt(variances(:,i-1)*delta_t).*rand2;
end

discount  = exp(-rate*option_maturity);
payoffs =  discount * max((stock_prices(:,num_periods+1)-strike),0);
Ms= [100,500,1000,5000,10000,30000,50000,70000,100000];
payoff_list = [mean(payoffs(1:100,:)),mean(payoffs(1:500,:)),mean(payoffs(1:1000,:)),mean(payoffs(1:5000,:)),mean(payoffs(1:10000,:)),mean(payoffs(1:30000,:)),mean(payoffs(1:50000,:)),mean(payoffs(1:70000,:)),mean(payoffs(1:100000,:))];

[call_price,call_delta] = bs_call(last_price,strike,rate,option_maturity,long_term_sigma,1);

plot(Ms,payoff_list,Ms,call_price*ones(1,length(Ms)))
legend({'Monte Carlo Simulation with Stochastic Vol','Black-Scholes Price'})

BSCallPrices = zeros(M,num_periods+1);
BSCallPrices(:,1) = call_price;

Deltas = zeros(M,num_periods+1);
Deltas(:,1) = call_delta;

PnL = zeros(M,num_periods+1);

for i=2:(num_periods+1)
    time_to_mat = option_maturity*((num_periods+1-i)/num_periods);
    passed_time = option_maturity - time_to_mat;
    [BSCallPrices(:,i),Deltas(:,i)] = bs_call(stock_prices(:,i),strike,rate,time_to_mat,sqrt(variances(:,i)),1);
    PnL(:,i) = PnL(:,i-1) + exp(rate*time_to_mat) * Deltas(:,i-1) .* (stock_prices(:,i) - stock_prices(:,i-1));
end

Final_PnL = -exp(int_rate*option_maturity) * BSCallPrices(:,1) + PnL(:,num_periods+1) + BSCallPrices(:,num_periods+1);
returns = Final_PnL./BSCallPrices(:,1);

mean(returns)
std(returns)
hist(returns)
