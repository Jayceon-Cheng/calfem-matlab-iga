clc;
clear;
close all;

file_name1 = '环形管.stp';

[lines_str1] = get_line(file_name1);

%将textB中的1到10行写到test_save文档中
for i = 222 : 355
    S = lines_str1{i};  
    fid = fopen('test_sX1.txt','at+');
        fprintf(fid,'%s\n',S);
     fclose(fid);
end
for i = 1635 : length(lines_str1)-1
    S = lines_str1{i};  
    fid = fopen('test_sX2.txt','at+');
        fprintf(fid,'%s\n',S);
     fclose(fid);
end
file_name2 ='test_sX1.txt';
[lines_str2] = get_line(file_name2);

file_name3 ='test_sX2.txt';
[lines_str3] = get_line(file_name3);
delete 'test_sX1.txt'  'test_sX2.txt'
p=0;
surf=struct('num',{},'order',{},'node',{});
for i = 1 :length( lines_str2)
    S = lines_str2{i};
    

    
        if S(1)=='#'
            p=p+1;s=0;
        pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      surf(p).num=regexp(S,pat,'match');
           
        elseif double(S(1))>=48&double(S(1))<=57
          pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';
      
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      surf(p).order=regexp(S,pat,'match');
        elseif double(S(1))==40
            s=s+1;
%          pat='[1-9][0-9]{0,3}';
           pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';   %浮点数+整数
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
        surf(p).node{s}=regexp(S,pat,'match');   
       
       end
       
end

Su=surf;


pp=0;
for i = 1 :length( lines_str3)
    S = lines_str3{i};
      
        if S(1)=='#'
            pp=pp+1;
%         pat='([0-9]*\.[0-9]+|[0-9]+)';
        pat='[-+]?[\d]+([\.][\d]*)?([Ee][+-]?[\d]+)|[-+]?([0-9]*\.[0-9]+|\d+)'; %浮点数+科学计数法
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      Node{pp,1}=regexp(S,pat,'match');
           
       end
       
end


%;
for i=1:3:length(surf)
   Ele=Su(i);
   for s=1:9
    Numb1(i,:,s)=str2num(char(Su(i).node{s})); 
   end
   p=0;
   for s=10:13
    p=p+1;
    knot1{i,:,p}=str2num(char(Su(i).node{s})); 
   end 
   
   m=0;
   for s=14:length(Su(i).node)
    m=m+1;
    weight1(i,:,m)=str2num(char(Su(i).node{s})); 
   end   
end
ii=0;
for i=2:3
    ii=ii+1;
   Ele=Su(i);
   for s=1:9
    Numb2(ii,:,s)=str2num(char(Su(i).node{s})); 
   end
   p=0;
   for s=10:13
    p=p+1;
    knot2{ii,:,p}=str2num(char(Su(i).node{s})); 
   end 
   
   m=0;
   for s=14:length(Su(i).node)
    m=m+1;
    weight2(ii,:,m)=str2num(char(Su(i).node{s})); 
   end   
end


% Numb= permute(Numb1,[1,3,2]);
B=cell(size(Numb2));
for i=1:numel(Numb2)
   a=Numb2(i);
   for s=1:size(a)
    for j=1:size(Node,1)
       lin= str2num(char(Node{j}));
       if a(s)==lin(1)
       B{i}(1:3)=lin(2:end);
       end
    end
   end
   B{i}(4)=weight2(i);
end
Order=[];

Order=[Order Su(1:4).order];
order=unique(str2num(char(Order)));

 Xi=[];
 Eta=[];
 Zeta=[];
 sum1=0;sum2=0;%size(knot1{1,1,3})
 for i=1:size(knot2{1,1,3})
   x1(1)=0;x2(1)=0;
   x1(i+1)=knot2{1,1,1}(i); x2(i+1)=knot2{1,1,2}(i); 
   sum1=sum1+x1(i+1);sum2=sum2+x2(i+1);
   Zeta((sum1-x1(i+1)+1):sum1)= knot2{1,1,3}(i);
   Eta((sum2-x2(i+1)+1):sum2)=knot2{1,1,4}(i);
 end
  sum3=0;
  for i=1:size(knot1{1,1,2},1)
   x3(1)=0;
   x3(i+1)=knot1{1,1,2}(i); 
   sum3=sum3+x3(i+1);
   Xi((sum3-x3(i+1)+1):sum3)=knot1{1,1,4}(i);
 end
 
 BJ=cell2mat(B);
%  B0=B;
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);

n = size(B,1);m = size(B,2);l = size(B,3);
% Order of basis
p = k_-n-1;q = l_-m-1;r = m_ -l -1;
KV0.Xi=Xi; KV0.Eta=Eta; KV0.Zeta=Zeta;


