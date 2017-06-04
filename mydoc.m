import mlreportgen.dom.*
d = Document('mydoc','docx');



%% insert image 
f=['figure_5','.jpg'];
imageObj = Image(f);
imageObj.Style={ScaleToFit};
append(d,imageObj)

%% figure 2
a=['figure_2','.jpg'];
imageObj = Image(a);
imageObj.Style={ScaleToFit};
append(d,imageObj)

%% insert image 2
b=['figure_1','.jpg'];
imageObj = Image(b);
imageObj.Style={ScaleToFit};
append(d,imageObj)

%% insert image 3
c=['figure_3','.jpg'];
imageObj = Image(c);
imageObj.Style={ScaleToFit};
append(d,imageObj)

%% insert image 4 

e=['figure_4','.jpg'];
imageObj = Image(e);
imageObj.Style={ScaleToFit};
append(d,imageObj)

%% image 5
close(d);
rptview('mydoc','docx');







% a = magic(5);
% [v,i] = max(a);
% [v1,i1] = max(max(a));
% table = Table(a);
% 
% text = table.entry(i(i1),i1).Children(1);
% text.Color = 'red';
% append(d,table);
% rank = 5;
% table = Table(magic(rank));
% table.Border = 'single';
% table.BorderWidth = '1px';
% 
% grps(1) = TableColSpecGroup;
% grps(1).Span = rank;
% grps(1).Style = {Width('0.2in'),Color('green')};
% 
% specs(1) = TableColSpec;
% specs(1).Span = 1;
% specs(1).Style = {Width('0.5in'),Bold,Color('red')};
% 
% grps(1).ColSpecs = specs;
% 
% table.ColSpecGroups = grps;
% append(d,table);

% %% write table 
% 
% A = [mu(1,1);sigma(1,1);p(1,1);q(1,1)];
% B = [mu(1,2);sigma(1,2);q(1,2);p(1,2)];
% P = [moyenne(1,1);ecart_type(1,1);qgauche(1,1);qdroite(1,1)];
% Q = [moyenne(1,2);ecart_type(1,2);qgauche(1,2);qdroite(1,2)];
% 
% Nom = {'A';'B';'P';'Q'};
% t = FormalTable({Nom,A,B,P,Q}); 
% 
% 
% r = TableRow();
% append(r,TableEntry('Paramètres'))
% append(r,TableEntry('Loi1final'));
% append(r,TableEntry('Loi2final'));
% append(r,TableEntry('Loi1ini'));
% append(r,TableEntry('Loi2ini'));
% append(t.Header,r);
% t.Style = {Border('inset','black','3px'), ...
%                ColSep('single','black','1px'), ...
%                RowSep('single','black','3px')};
% append(d,t);

% %% write table 2
% 
% Loi_1 = [pbcriseL1;pbnocriseL1];
% Loi_2 = [pbcriseL2;pbnocriseL2];
% 
% Nom = {'Loi_1';'Loi_2'};
% tt = FormalTable({Nom,Loi_1,Loi_2}); 
% rr = TableRow();
% append(rr,TableEntry('Etat'))
% append(rr,TableEntry('Crise'));
% append(rr,TableEntry('Non-Crise'));
% append(t.Header,rr);
% tt.Style = {Border('inset','black','3px'), ...
%                ColSep('single','black','1px'), ...
%                RowSep('single','black','3px')};
% append(d,tt);