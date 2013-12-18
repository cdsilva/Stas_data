clear all
close all

res = '-r300';
fmt = '-djpeg';

set(0,'DefaultAxesFontSize',14)

%% toy example 1
theta = 0:0.1:pi;
theta2 = [1 pi/2 pi-1];
amp = 0.1:0.1:1;

m = length(amp);
n = length(theta);
n2 = length(theta2);

data = (amp' * ones(1,n)) .* (ones(m,1) * sin(theta).^4);
data2 = (amp' * ones(1,n2)) .* (ones(m,1) * sin(theta2).^4);
time = 1:m;

figure;
imagesc(data)
xlabel('position')
ylabel('time')
print('example1_all',fmt,res)

figure;
imagesc(data2)
xlabel('position')
ylabel('time')
print('example1_coarse',fmt,res)

figure;
hold on
for i=1:length(amp)
    plot(theta, amp(i)*sin(theta).^4)
    plot(theta2, amp(i)*sin(theta2).^4,'o')
end

figure;
scatter3(data2(:,1),data2(:,2),data2(:,3),500,time,'.')
xlabel('position 1')
ylabel('position2')
zlabel('position3')
grid on
print('example1_embedding',fmt,res)

%% dmaps example
s=0.1:0.1:5;

a=2;
b=0.2;

t=(1/b)*log(b*s/(a*sqrt(1+b^2)));

x=a*exp(b*t).*cos(t);
y=a*exp(b*t).*sin(t);

X=[x; y; zeros(1,length(t))];
[N,n]=size(X);

%rotate first data set  by angle 1
angle=pi/2-0.1;

R=[1 0 0;
    0 cos(angle) -sin(angle);
    0 sin(angle) cos(angle)];

X1=R*X;

W = squareform(pdist(X1'));
W = W.^2;
[V, D] = dmaps(W, 0.05, 10);

figure;
scatter3(X1(1,:),X1(2,:),X1(3,:),500,V(:,2),'.')
xlabel('position 1')
ylabel('position2')
zlabel('position3')
grid on
print('dmaps_spiral',fmt,res)
