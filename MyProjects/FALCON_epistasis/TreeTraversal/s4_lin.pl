#!/usr/bin/perl -w
open(INP, $ARGV[0]); 
open(OUT, ">xx.txt"); 


@a=<INP>;
close INP;

$xx=0;
$xxx=0;
$zz=0;
$zzz=0;
$rrx=0;
$rxx=0;
for ($u=0; $u<@a; $u+=50){

@p=();
for ($t=0; $t<5; $t++){
for ($tt=0; $tt<5; $tt++){
for ($ttt=0; $ttt<5; $ttt++){
$p[$t][$tt][$ttt]=0;
}
}
}

$p[0][0][0]=1;

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[1]==1)&&($b[0]==0)){
$p[1][0][0]=$c[1];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[1]==2)&&($b[0]==0)){
$p[2][0][0]=$c[1];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[1]==3)&&($b[0]==0)){
$p[3][0][0]=$c[1];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[1]==4)&&($b[0]==0)){
$p[4][0][0]=$c[1];
last;
}
}












for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==2)&&($b[1]==1)&&($b[0]==0)){
$p[1][2][0]=$c[2];
$p[2][1][0]=$c[2];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==3)&&($b[1]==1)&&($b[0]==0)){
$p[1][3][0]=$c[2];
$p[3][1][0]=$c[2];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==3)&&($b[1]==2)&&($b[0]==0)){
$p[2][3][0]=$c[2];
$p[3][2][0]=$c[2];
last;
}
}


for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==4)&&($b[1]==1)&&($b[0]==0)){
$p[1][4][0]=$c[2];
$p[4][1][0]=$c[2];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==4)&&($b[1]==2)&&($b[0]==0)){
$p[2][4][0]=$c[2];
$p[4][2][0]=$c[2];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[2]==4)&&($b[1]==3)&&($b[0]==0)){
$p[3][4][0]=$c[2];
$p[4][3][0]=$c[2];
last;
}
}























for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[3]==3)&&($b[2]==2)&&($b[1]==1)&&($b[0]==0)){
$p[1][2][3]=$c[3];
$p[1][3][2]=$c[3];
$p[3][2][1]=$c[3];
$p[3][1][2]=$c[3];
$p[2][1][3]=$c[3];
$p[2][3][1]=$c[3];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[3]==4)&&($b[2]==2)&&($b[1]==1)&&($b[0]==0)){
$p[1][2][4]=$c[3];
$p[1][4][2]=$c[3];
$p[4][2][1]=$c[3];
$p[4][1][2]=$c[3];
$p[2][1][4]=$c[3];
$p[2][4][1]=$c[3];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[3]==4)&&($b[2]==3)&&($b[1]==1)&&($b[0]==0)){
$p[1][3][4]=$c[3];
$p[1][4][3]=$c[3];
$p[4][3][1]=$c[3];
$p[4][1][3]=$c[3];
$p[3][1][4]=$c[3];
$p[3][4][1]=$c[3];
last;
}
}

for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[3]==4)&&($b[2]==3)&&($b[1]==2)&&($b[0]==0)){
$p[2][3][4]=$c[3];
$p[2][4][3]=$c[3];
$p[4][3][2]=$c[3];
$p[4][2][3]=$c[3];
$p[3][2][4]=$c[3];
$p[3][4][2]=$c[3];
last;
}
}



for ($v=2; $v<50; $v+=2){
chomp $a[$u+$v];
@b=();
@b=split(/	/, $a[$u+$v]);
chomp $a[$u+$v+1];
@c=();
@c=split(/	/, $a[$u+$v+1]);
if (($b[4]==4)&&($b[3]==3)&&($b[2]==2)&&($b[1]==1)&&($b[0]==0)){
$pp=$c[4];
last;
}
}









for ($j=1; $j<5; $j++){

@mean=();
@slope=();

$x=0;


for ($k=1; $k<5; $k++){
if ($k!=$j){
@d=();
@e=();

$e[0]=(($p[$j][0][0]-$p[0][0][0])/$p[0][0][0]);
$d[0]=$p[0][0][0];
$e[1]=(($p[$j][$k][0]-$p[$k][0][0])/$p[$k][0][0]);
$d[1]=$p[$k][0][0];

$mean[$x]=$p[$k][0][0];


$jd=0;
$je=0;
$kc=0;

for ($tt=0; $tt<@d; $tt++){

chomp $d[$tt];
$jd=($jd+$d[$tt]);

chomp $e[$tt];
$je=($je+$e[$tt]);

$kc++;
}


$pd=($jd/$kc);
$pe=($je/$kc);


$ah=0;
$bh=0;
for ($hh=0; $hh<@d; $hh++){
$ah+=(($d[$hh]-$pd)*($e[$hh]-$pe));
$bh+=(($d[$hh]-$pd)*($d[$hh]-$pd));
}
if ($bh!=0){
$slope[$x]=($ah/$bh);
$x++;
}


}
}




@f=('1','2','3','4');
@f=grep(/[^$j]/, @f);

for ($x=0; ($x+1)<@f; $x++){
$k=$f[$x];
for ($y=($x+1); $y<@f; $y++){
$h=$f[$y];

@d=();
@e=();
if ($p[$k][$h][0]==0){
}
else {
$e[0]=(($p[$j][0][0]-$p[0][0][0])/$p[0][0][0]);
$d[0]=$p[0][0][0];
$e[1]=(($p[$j][$k][0]-$p[$k][0][0])/$p[$k][0][0]);
$d[1]=$p[$k][0][0];
$e[2]=(($p[$j][$h][0]-$p[$h][0][0])/$p[$h][0][0]);
$d[2]=$p[$h][0][0];
$e[3]=(($p[$j][$k][$h]-$p[$k][$h][0])/$p[$k][$h][0]);
$d[3]=$p[$k][$h][0];


$mean[$x]=(($p[$k][0][0]+$p[$h][0][0])/2);


$jd=0;
$je=0;
$kc=0;

for ($tt=0; $tt<@d; $tt++){

chomp $d[$tt];
$jd=($jd+$d[$tt]);

chomp $e[$tt];
$je=($je+$e[$tt]);

$kc++;
}


$pd=($jd/$kc);
$pe=($je/$kc);


$ah=0;
$bh=0;
for ($hh=0; $hh<@d; $hh++){
$ah+=(($d[$hh]-$pd)*($e[$hh]-$pe));
$bh+=(($d[$hh]-$pd)*($d[$hh]-$pd));
}
if ($bh!=0){
$slope[$x]=($ah/$bh);
$x++;
}

}
}
}










@g=('1','2','3','4');
@g=grep(/[^$j]/, @g);
$k=$g[0];
$h=$g[1];
$l=$g[2];

@d=();
@e=();

if (($p[$k][$h][$l]==0)||($p[$h][$l][0]==0)||($p[$k][$l][0]==0)||($p[$k][$h][0]==0)){
}
else {
$e[0]=(($p[$j][0][0]-$p[0][0][0])/$p[0][0][0]);
$d[0]=$p[0][0][0];
$e[1]=(($p[$j][$k][0]-$p[$k][0][0])/$p[$k][0][0]);
$d[1]=$p[$k][0][0];
$e[2]=(($p[$j][$h][0]-$p[$h][0][0])/$p[$h][0][0]);
$d[2]=$p[$h][0][0];
$e[3]=(($p[$j][$l][0]-$p[$l][0][0])/$p[$l][0][0]);
$d[3]=$p[$l][0][0];
$e[4]=(($p[$j][$k][$h]-$p[$k][$h][0])/$p[$k][$h][0]);
$d[4]=$p[$k][$h][0];
$e[5]=(($p[$j][$k][$l]-$p[$k][$l][0])/$p[$k][$l][0]);
$d[5]=$p[$k][$l][0];
$e[6]=(($p[$j][$h][$l]-$p[$h][$l][0])/$p[$h][$l][0]);
$d[6]=$p[$h][$l][0];
$e[7]=(($pp-$p[$k][$h][$l])/$p[$k][$h][$l]);
$d[7]=$p[$k][$h][$l];

$mean[$x]=(($p[$k][0][0]+$p[$h][0][0]+$p[$l][0][0])/3);


$jd=0;
$je=0;
$kc=0;

for ($tt=0; $tt<@d; $tt++){

chomp $d[$tt];
$jd=($jd+$d[$tt]);

chomp $e[$tt];
$je=($je+$e[$tt]);

$kc++;
}


$pd=($jd/$kc);
$pe=($je/$kc);


$ah=0;
$bh=0;
for ($hh=0; $hh<@d; $hh++){
$ah+=(($d[$hh]-$pd)*($e[$hh]-$pe));
$bh+=(($d[$hh]-$pd)*($d[$hh]-$pd));
}
if ($bh!=0){
$slope[$x]=($ah/$bh);
$x++;
}


}







$jx=0;
$jy=0;
$kx=0;

for ($tx=0; $tx<@mean; $tx++){

chomp $mean[$tx];
$jx=($jx+$mean[$tx]);

chomp $slope[$tx];
$jy=($jy+$slope[$tx]);

$kx++;
}


$px=($jx/$kx);
$py=($jy/$kx);


$ax=0;
$bx=0;
$cx=0;
for ($hx=0; $hx<@mean; $hx++){
$ax+=(($mean[$hx]-$px)*($slope[$hx]-$py));
$bx+=(($mean[$hx]-$px)*($mean[$hx]-$px));
$cx+=(($slope[$hx]-$py)*($slope[$hx]-$py));
}
if (($bx*$cx)!=0){
$rx=($ax/sqrt($bx*$cx));
if ($rx<0){

$zz++;
@sort1=();
@sort1=sort { $a <=> $b } @mean;
for ($bc=0; $bc<@mean; $bc++){
if ($mean[$bc]==$sort1[0]){
if ($slope[$bc]>0){
$zzz++;
}
}
}

$xx++;
$rrx+=$rx;
}
if ($rx>=0){
$rxx+=$rx;
$xxx++;
}
}


}

}

$yyy=($xx/($xx+$xxx));
$uu=($rrx/$xx);
$uuu=(($rxx+$rrx)/($xx+$xxx));
$yyyy=($zzz/$zz);
print OUT "Percentage of negative correlation in s1.txt output is $yyy", "\n";
print OUT "mean of all negative correlation coefficients in s1.txt output is $uu", "\n";
print OUT "mean of all correlation coefficients in s1.txt output is $uuu", "\n";
print OUT "Percentage of positive AP as the first data point in s3.txt output is $yyyy", "\n";
close OUT;
