close all
clear all

A=    1.0e-09 *[ 0.1818
    0.2168
    0.1951
    0.1988
    0.2019
    0.1943
    0.2029
    0.2017
    0.1976
    0.1917
    0.2008
    0.2119
    0.2129
    0.1961
    0.1939
    0.2080
    0.1938
    0.2141
    0.2035
    0.1969
    0.1970
    0.2151
    0.2024
    0.1986
    0.2131
    0.1883
    0.2090
    0.2084
    0.2011
    0.1924
    0.2016
    0.2040
    0.1980
    0.1992
    0.1952
    0.1959
    0.2050
    0.1943
    0.2063
    0.2001
    0.2136
    0.1886
    0.1915
    0.1829
    0.2070
    0.2031
    0.2140
    0.2040
    0.2030
    0.1893
    0.2217
    0.2022
    0.1881
    0.2048
    0.1932
    0.2062
    0.1984
    0.2206
    0.1921
    0.2182
    0.1945
    0.1945
    0.2046
    0.1917
    0.1940
    0.2030
    0.2038
    0.2105
    0.2084
    0.1831
    0.1996
    0.1939
    0.2157
    0.1958
    0.2080
    0.2006
    0.2133
    0.1877
    0.2037
    0.1871
    0.1859
    0.1962
    0.2212
    0.1916
    0.1910
    0.2008
    0.2033
    0.2041
    0.2058
    0.2138
    0.2011
    0.1937
    0.2172
    0.1929
    0.2189
    0.2008
    0.2035
    0.1942
    0.2011
    0.1894];

B= [
    2.3000
    2.6303
    2.9606
    3.2909
    3.6212
    3.9515
    4.2818
    4.6121
    4.9424
    5.2727
    5.6030
    5.9333
    6.2636
    6.5939
    6.9242
    7.2545
    7.5848
    7.9152
    8.2455
    8.5758
    8.9061
    9.2364
    9.5667
    9.8970
   10.2273
   10.5576
   10.8879
   11.2182
   11.5485
   11.8788
   12.2091
   12.5394
   12.8697
   13.2000
   13.5303
   13.8606
   14.1909
   14.5212
   14.8515
   15.1818
   15.5121
   15.8424
   16.1727
   16.5030
   16.8333
   17.1636
   17.4939
   17.8242
   18.1545
   18.4848
   18.8152
   19.1455
   19.4758
   19.8061
   20.1364
   20.4667
   20.7970
   21.1273
   21.4576
   21.7879
   22.1182
   22.4485
   22.7788
   23.1091
   23.4394
   23.7697
   24.1000
   24.4303
   24.7606
   25.0909
   25.4212
   25.7515
   26.0818
   26.4121
   26.7424
   27.0727
   27.4030
   27.7333
   28.0636
   28.3939
   28.7242
   29.0545
   29.3848
   29.7152
   30.0455
   30.3758
   30.7061
   31.0364
   31.3667
   31.6970
   32.0273
   32.3576
   32.6879
   33.0182
   33.3485
   33.6788
   34.0091
   34.3394
   34.6697
   35.0000];
plot(B,A,'*', 'Linewidth', 4)
hold on
xlim([0 35])
idx = find(min(abs(B-23))==abs(B-23));
plot(B(idx),A(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
set(gca, 'XTick', [0:2:35])
xlabel('$g\_Na$','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_Na)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
shg