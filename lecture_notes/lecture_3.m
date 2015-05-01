%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test MATLAB version and installed toolboxes
ver

v = ver('matlab');
v.Version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% warming up: matrix manipulations

A = [1 2 4 8; 8 16 32 64; 64 128 256 512]

sum(A)

%% sum in columns - general Matlab feature!

sum(A')
sum(A,2)

sum(sum(A))
sum(A(:))

A = [-1 -3 -5; 1 3 4]
B = [2 3 1; -4 1 4]

%% concatenation
cat(1,A,B)
cat(2,A,B)

vertcat(A,B)
horzcat(A,B)

C = cat(3,A,B)

size(C)

%% flipping
A = [1 2 4 8; 8 16 32 64; 64 128 256 512]
fliplr(A)
flipud(A)

permute(A,[2 1])
permute(A,[3 2 1])
A

squeeze(A)

%% circularly shifts elements in matrix
v = ver('matlab');
if str2num(v.Version) < 8.4
   circshift(A,1)
   circshift(A',1)'

else
% new matlab
   circshift(A,1,2)
end


%% reshape
B = reshape(A,1,numel(A))
B = reshape(A,numel(A),1)

C = reshape(A,2,6)
C = reshape(A,6,2)

%% matrix elements/regions
A = rand(5,5)
tril(A)
triu(A)

diag(A)

eye(3,3)



%% sorting
A = rand(3,6)

sort(A)
sort(A,2)
sort(A,2,'descend')

sortrows(A,2)

%% comparison
A = [-1 2 -5; 1 6 4]
B = [2 3 1; -4 1 4]

ismember(A,B)

unique(B)
union(A,B)
intersect(A,B)
A(ismember(A,B))
setdiff(A,B)
setdiff(B,A)


A > 0
A(A>0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sparse matrices
% sparse matrices are an additional data class
% very useful when (spatial) data contains many zeros
A = [-1 2 0 -5; 0 9 2 0; 8 0 6 4]
B = sparse(A)

whos

speye(3,3)
spones(B)
C = full(B)

X = [1 1 0.6; 1 2 0.9; 1 3 0.3; 
     2 1 0.2; 2 2 0.8; 2 3 0.4]
spconvert(X)

sprand(B)
sprand(3,3,.4)
spy(sprand(10,10,.5))
spy(sprandsym(10,.5))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preallocating memory, vectorisation
% 

%% 1st example without memory preallocation
% the matrix x is growing with each iteration step and, thus,
% requires new, time consuming memory preallocation with each iteration step
clear
N = 2500;
tic
for i = 1:N
   for j = 1:N
      x(i,j) = i^2 * j;
   end
end
toc


%% 2nd example with memory preallocation
clear
N = 2500;
tic
x = zeros(N,N);
for i = 1:N
   for j = 1:N
      x(i,j) = i^2 * j;
   end
end
toc

%% vectorisation 1
% MATLAB allows vector operations; can help to avoid
% time consuming loops.
% However, it is important in which order we vectorise!
% In the following example, the vectorisation is over the
% 2nd dimension (row), but takes long time.
clear
N = 5500;
tic
x = zeros(N,N);
j = 1:N;
for i = 1:N
   x(i,:) = i^2 * j;
end
toc

%% vectorisation 2
% In the following example, the vectorisation is over the
% 1st dimension (colum), causing significant speedup!
clear
N = 5500;
tic
x = zeros(N,N);
j = (1:N)';
for i = 1:N
   x(:,i) = i^2 * j;
end
toc

%% complete vectorisation
% even complete vectorisation of both loops is possible,
% can help to safe further time, but not in every case
% (depends on the operations and the matrix size)
clear
N = 5500;
tic
x = zeros(N,N);
i = (1:N)';
j = 1:N;
x = repmat(i,1,N).^2 .* repmat(j,N,1);
toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% floating point issues
% precission
realmin('double')
realmax('double')

eps(1)

sin(pi)

format long
format hex
format bank
format


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I/O
%% simple import
x = load('example.dat');
plot(x(:,2),x(:,3))

%%
% does not work in every case
x = load('tsi_DB_lin_40_11yr.txt')

%% using importdata
% importdata is a very intelligent importing function, working in
% most cases
x = importdata('tsi_DB_lin_40_11yr.txt')

%% 
% variable x is now a structure
plot(x.data(:,1),x.data(:,2))
xlabel(x.colheaders{1})
ylabel(x.colheaders{2})

%% using textscan
% quite flexible but rather complex settings
fid = fopen('tsi_DB_lin_40_11yr.txt');
data = textscan(fid,'%u %f %f','Headerlines',4,'collectoutput',1)
fclose(fid);

year = data{1};
DB = data{2};
plot(year,DB(:,1))

%% MATLAB date functions
% MATLAB offers a convenient way to handle date information
x = importdata('Mawmluh_village_precp2011.dat')
x.textdata(:,1)
time = x.textdata(2:end,1);

%%
% datenum converts a date into a number that counts the number of
% days after year 0000
datenum(time{1},'dd.mm.yy')

%%
% datestr produces formatted output from the datenum value
datestr(datenum(time{1},'dd.mm.yy'))

t = datenum(time,'dd.mm.yy')

bar(t,x.data)

%%
% datetick changes the axis labels to formatted strings
datetick

datetick('x',2)
datetick('x','mm')

%% data import using fgetl
% for more complex files, where importdata and textscan will not work
clear all
fid = fopen('Mawmluh_village_precp2011.dat');
cnt = 0;
while 1 % infinite loop
   temp = fgetl(fid); % reads line per line
   if ~ischar(temp) % infinite loop will stop when end of file is reached (indicated when TEMP contains a number instead of a string)
      break
   end
   if cnt > 0
      [txt1 rem1] = strtok(temp,';'); % allows separating strings by a given delimiter, here ";"
      t(cnt) = datenum(txt1,'dd.mm.yy');
      data(cnt) = str2num(rem1);
   end
   cnt = cnt + 1;
end
fclose(fid);

bar(t,data)
datetick

%% using fread
% intersting for binary data or extremely complex data;
% finally needs intensive coding for parsing the information
fid = fopen('Mawmluh_village_precp2011.dat');
temp = fread(fid);
fclose(fid);

temp
char(temp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exporting data

x = load('example.dat');

% convert year in BP!
data = [1950-x(:,2) x(:,3)];

%% export using save
save example2.dat data -ascii -tabs

%% export using fopen and fprintf
% fprintf allows nicely formatted output
fid = fopen('example2.dat','w');
fprintf(fid,'Year (AD)\td18O\n');
fprintf(fid,'%i\t%6.4f\n',data');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sliding windows and convolution
clear
x = importdata('Mawmluh_village_precp2011.dat')
data = x.data;
t = datenum(x.textdata(2:end,1),'dd.mm.yy')

%% 14 days moving average (w/o overlap)
mData = [];
cnt = 1;
for i = 1:14:length(t)-14
   mData(cnt) = mean(data(i:i+13));
   cnt = cnt + 1;
end
plot(t,data,t(1:14:end-14),mData)
datetick

%% 14 days moving average (w/ overlap)
mData = [];
for i = 1:length(t)-14
   mData(i) = mean(data(i:i+13));
end
plot(t,data,t(1:end-14),mData)
datetick

%% convolution with rectangular filter = moving average
b = ones(14,1)/14;

mData = conv(data,b);
whos

%%
% ensuring the same length of the resulting vector
mData = conv(data,b,'same');
whos

plot(t,data,t,mData)
datetick


%% convolution with Gaussian kernel filter 
b = gausswin(14);
plot(b)

mData = conv(data,b,'same');

plot(t,data,t,mData)
datetick


%% windowed analysis using buffer
% w/o overlap
wData = buffer(data,14);
mData = mean(wData)

plot(t,data,t(1:14:end),mData)
datetick

%%
% w/ overlap (using 7 days overlap)
wData = buffer(data,14,7);
mData = mean(wData)

plot(t,data,t(1:7:end),mData)
datetick



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parallelisation
% required: multicore CPU, Parallel Computing Toolbox
n = 10000;
x = randn(1,n);
y = zeros(1,n);
matlabpool('open','12');

tic
for i = 1 : n
    y(i) = std(x(1:i));
end
t1 = toc;
tic
parfor i = 1 : n
    y(i) = std(x(1:i));
end
t2 = toc;

fprintf('\n Normal for: %f secs',t1);
fprintf('\n     parFor: %f secs\n\n',t2);
matlabpool('close');

%%
% Result on Xeon workstation with 16 cores:
% Normal: 1.58 sec, Parfor: 0.33 sec
% 
% Remark: several builtin function already use multi-threading
% automatically, parfor therefore not as powerful for such
% functions, but useful for own functionss


%% computation using GPU
% the latest versions of Parallel Computing Toolbox allow computations 
% on the GPU

%%
% example from beginning of the lecture
clear
N = 5500;
tic
x = zeros(N,N);
j = (1:N)';
for i = 1:N
   x(:,i) = i^2 * j;
end
toc

%% 
% the same calculation on a GPU
clear
N = 5500;
tic
i = gpuArray((1:N)');
j = gpuArray(1:N);
x = repmat(i,1,N).^2 .* repmat(j,N,1);
toc

%%
% on a Xeon with Nvidia Quadro 5000
% normal: 0.18 sec, GPU: 0.05 sec


%% profiler
% allows for detailed reporting of bottlenecks in functions and scripts
profile on
ar2
profreport


