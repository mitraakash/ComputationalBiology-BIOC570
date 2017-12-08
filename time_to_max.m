function [time_course time] = time_to_max(a, x0)
N = 990;
rhs = @(t,x) [a*x(1)*(1-x(1))];
x = x0;
X = x/N;
sol1 = ode23(rhs, [x N], X);
time_course = (sol1.x);
time = length(time_course);
end
