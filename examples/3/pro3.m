
% clc;
% clear;
close all;

file_name1 = '近似螺丝刀2.stp';

[lines_str1] = get_line(file_name1);

%将textB中的1到10行写到test_save文档中
for i = 378 :504
    S = lines_str1{i};  
    fid = fopen('test_1.txt','at+');
        fprintf(fid,'%s\n',S);
    fclose(fid);
end
for i = 777 : length(lines_str1)-1
    S = lines_str1{i};  
    fid = fopen('test_2.txt','at+');
        fprintf(fid,'%s\n',S);
     fclose(fid);
end
file_name2 ='test_1.txt';
[lines_str2] = get_line(file_name2);

file_name3 ='test_2.txt';
[lines_str3] = get_line(file_name3);
delete 'test_1.txt' 'test_2.txt'
p=0;nn=0;pj=0;ni=1;
surf=struct('num',{},'order',{},'con_Node',{});

%------提取所需输入的信息并进行初步归类---------
for i = 1 :length( lines_str2)
 
    S = lines_str2{i};
    if (i>=2)
    Si = lines_str2{i-1};  
    end
        if S(1)=='#'
            p=p+1;ni=1;s=0;
        pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      surf(p).num=regexp(S,pat,'match');
           
        else
            if double(S(1))>=48&double(S(1))<=57
          pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';
      
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      surf(p).order=regexp(S,pat,'match');
            else
                if double(S(1))==40    
                  if double(Si(1))>=65&double(Si(1))<=90
                   ni=ni+1; s=0;
                  end
                  s=s+1;
%          pat='[1-9][0-9]{0,3}';
             pat='[-+]?([0-9]*\.[0-9]+|[0-9]+)';   %浮点数+整数
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
           surf(p).con_Node(ni).node{s}=regexp(S,pat,'match');  
                end 
            end        
        end       
end

%------提取每个点的坐标------
Su=surf;

pp=0;

for i = 1 :length( lines_str3)
    S = lines_str3{i};
      
        if S(1)=='#'
            pp=pp+1;
%         pat='([0-9]*\.[0-9]+|[0-9]+)';
        pat='[-+]?[\d]+([\.][\d]*)?([Ee][+-]?[\d]+)|[-+]?([0-9]*\.[0-9]+|\d+)'; %浮点数+科学计数法
%     pat='(([0-9]|#)?)\..*?(?(1)([1-9][0-9]{0,3}))\s'
      Coord{pp,1}=regexp(S,pat,'match');
           
       end
       
end

%{
for pi=1:p
 for i=1:length(surf(pi).con_Node(1).node)
   Numb(i,:,pi)=str2num(char(Su(pi).con_Node(1).node{i}));
 end
 for i=1:length(surf(pi).con_Node(2).node)
 Knot{i,:,pi}=str2num(char(Su(pi).con_Node(2).node{i}));
 end
 for i=1:length(surf(pi).con_Node(3).node)
 Weight(i,:,pi)=str2num(char(Su(pi).con_Node(3).node{i}));
 end
end
%}

%-----将提取的初步信息转换成可用的规范的形式-----
zn=0;
for i=(length(Su)/2+1):length(Su)    %提取模型的内外曲面的控制点
    Ele=Su(i); zn=zn+1;
   for j=1:length(surf(i).con_Node(1).node)
   Numb(:,j,zn)=str2num(char(Su(i).con_Node(1).node{j}));
   end
   
   for j=1:length(surf(i).con_Node(3).node)
   Weight(:,j,zn)=str2num(char(Su(i).con_Node(3).node{j}));
   end
end   
zn=0;
for i=1:2:length(Su)        %提取不同方向上的两个曲面的结向量情况
    Ele=Su(i); zn=zn+1;
   for j=1:length(surf(i).con_Node(2).node)
    Knot{:,j,zn}=str2num(char(Su(i).con_Node(2).node{j}));
   end
end

%{
for i=1:2
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
for i=3:4
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
%}

%-------将控制点信息包括坐标和权进行整合成B矩阵-----
B=cell(size(Numb));
for i=1:numel(Numb)
   a=Numb(i);
   for s=1:size(a)
    for j=1:size(Coord,1)
       lin= str2num(char(Coord{j}));
       if a(s)==lin(1)
       B{i}(1:3)=lin(2:end);
       end
    end
   end
   B{i}(4)=Weight(i);
end

%--------形函数的度数degree-----
Order=[];
Order=[Order Su.order];
order=unique(str2num(char(Order)));

% order=[order; 1];
%{
Xi=[];
 Eta=[];
 Zeta=[];
 sum1=0;sum2=0;%size(knot1{1,1,3})
 for i=1:size(knot1{1,1,3})
   x1(1)=0;x2(1)=0;
   x1(i+1)=knot1{1,1,1}(i); x2(i+1)=knot1{1,1,2}(i); 
   sum1=sum1+x1(i+1);sum2=sum2+x2(i+1);
   Zeta((sum1-x1(i+1)+1):sum1)= knot1{1,1,3}(i);
   Eta((sum2-x2(i+1)+1):sum2)=knot1{1,1,4}(i);
 end
  sum3=0;
  for i=1:size(knot2{1,1,2},1)
   x3(1)=0;
   x3(i+1)=knot2{1,1,2}(i); 
   sum3=sum3+x3(i+1);
   Xi((sum3-x3(i+1)+1):sum3)=knot2{1,1,4}(i);
  end
%}

%--------分配结向量--------
 Xi=[];
 Eta=[];
 Zeta=[];
 numknot=size(Knot) %[1,4,]
 sum1=0;sum2=0;%size(knot1{1,1,3})
 for i=1:numknot(3)
     for j=1:numknot(2)/2
         sum1=0;
      for k=1:length(Knot{:,j,i})
         x1(1)=0;
         x1(i+1)=Knot{:,j,i}(k);  
         sum1=sum1+x1(i+1);
         kv{i,j}((sum1-x1(i+1)+1):sum1)= Knot{:,j+2,i}(k);
      end
     end
 end
 
 zn=0;
for i=1:2:length(Su)        %提取不同方向上的两个曲面的结向量情况
    Ele=Su(i); zn=zn+1;
    deKV(zn,:)=str2num(char(Su(i).order))';
end

N=size(B);
for i=1:numel(kv)
     Num=length(kv{i})-deKV(i)-1;
     if Num ==N(1)
         Xi=kv{i}; p=deKV(i);   
     elseif Num ==N(2)
         Eta=kv{i}; q=deKV(i);  
     elseif Num ==N(3)
         Zeta=kv{i}; r=deKV(i);
     end
end

 %{
  if length(size(B))<=2
     Zeta=[1,1];
     for i=1:size(Knot{1,1,s})
     x1(1)=0;
     x1(i+1)=Knot{1,1,s}(i);  
     sum1=sum1+x1(i+1);
     Xi((sum1-x1(i+1)+1):sum1)= Knot{3,1,s}(i);
     end
     for i=1:size(Knot{2,1,s})
     x2(1)=0;
     x2(i+1)=Knot{2,1,s}(i); 
     sum2=sum2+x2(i+1);
      Eta((sum2-x2(i+1)+1):sum2)=Knot{4,1,s}(i);
     end
  else
   sum3=0;
   for i=1:length(knot2{1,1,1})
   x3(1)=0;
   x3(i+1)=knot2{1,1,2}(i); 
   sum3=sum3+x3(i+1);
   Zeta((sum3-x3(i+1)+1):sum3)=knot2{1,1,4}(i);
   end
   for i=1:size(Knot{:,1,s})
   x1(1)=0;
   x1(i+1)=Knot{1,1,s}(i);  
   sum1=sum1+x1(i+1);
   Xi((sum1-x1(i+1)+1):sum1)= Knot{3,1,s}(i);
  end
  for i=1:size(Knot{2,1,s})
  x2(1)=0;
  x2(i+1)=Knot{2,1,s}(i); 
   sum2=sum2+x2(i+1);
   Eta((sum2-x2(i+1)+1):sum2)=Knot{4,1,s}(i);
  end
  end
%}
%  
% B=B(:,2:size(B,2),1);
 BJ=cell2mat(B);
 KV0.Xi=Xi; KV0.Eta=Eta; KV0.Zeta=Zeta;
 
% k_ = length(Xi);
% l_ = length(Eta);
% m_ = length(Zeta);
% 
% n = size(B,1);m = size(B,2);l = size(B,3);
% % Order of basis
% p = k_-n-1;q = l_-m-1;r = m_ -l -1;


