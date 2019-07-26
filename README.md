# Overfitting Least Squares
I stumbled upon this when I was doing a normal least squares fitting. I fitted the polynomial to the data but then asked myself: 
"What would the coefficients look like if I keep fitting this data to higher and higher degree polynomials." My initial thought was that the
sequence of coefficients for progressively larger polynomials
> a_n = (a_1, a_2, ..., a_n)

would converge to some true sequence which fits the data best. However, I was also trying to measure if that idea could be true with the inevitable Runge's Phenomenon:
## Runge's Phenomenon
Runge's Phenomenon is observed when overfitting and is described as a loss of smoothness. In the case of polynomials, erratic undulations form at the ends of the fit.

For example, this is what I would call a great fit to the data:

![perfectfit](https://user-images.githubusercontent.com/44763636/61954756-3210f880-afc2-11e9-8e6d-2b51e537dde7.png)

And this is overfitting with Runge's Phenomenon: 
> Here the first 100 polynomials fitted, from 
p1 = (a_1)x to 
p100 = (b_1)x+...+(b_100)x^100. 
Note how they converge to the same [overfit] expression

![overfit](https://user-images.githubusercontent.com/44763636/61954965-b794a880-afc2-11e9-94e3-1c0a15b1a074.png)

## So What?
So what does this mean? Well, despite Runge's Phenomenon, increasing polynomials will converge to some "best fit" (though practically not helpful) of the data. This led me to believe that the sequence of coefficients of each fit polynomial will converge. And, that is exactly what I found:

![1_n100](https://user-images.githubusercontent.com/44763636/61971327-7b277380-afe7-11e9-84c8-232fbf003191.png)
> The data used to create this is the same as the data used above

Let's take a look at that heat map of the coefficients:

![HeatMap](https://user-images.githubusercontent.com/44763636/61971310-6c40c100-afe7-11e9-9b79-05857107b9ae.png)
> Each *row* in the heatmap represents one polynomial, starting with the degree one polynomial on the top and ending with the degree 100 polynomial on the bottom. Each *column* represents the magnitude of the coefficient corresponding to that term in the polynomial. For example, the first column corresponds to the size of the coefficient on the degree 1 term of every polynomial. Consequently, this is why the top right half of the heatmap is all zero filled; a degree n polynomial has all coefficients on its terms of degree greater than n uniformly equal to 0.

This is interesting because it shows slow and inconsistent convergence.

Let's extract these motifs. To do so, I will use a motif finding algorithm (just for fun). So first we need to take the sequences of coefficients and turn them sequences of letters. To do this I will place them into 26 bins based on their magnitude, then assign each bin a letter of the alphabet. So for each coefficient in a sequence of coefficients, it is assigned a letter, this effectively changes the sequence of coefficients into a sequence of letters. Here is a snapshot of that for polynomials 20-29:

```
>20
nooopoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>21
nooppnooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>22
nooppnooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>23
nooppnooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>24
nonpqnnoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>25
nonpqnnpooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>26
nonpqnmppoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>27
nonpqomopoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>28
nonpqonopoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
>29
nonpqnnopoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\
```
> (This is fasta format)

Now I will run "meme" on it to get the motif. "meme" is part of the MEME Suite, a motif-based sequence analysis tool which can detect denovo motifs, which is specifically made for DNA or RNA sequences. However I can define a new alphabet and give it as an input.

Running **translateMotif.py** on the output gives the motif, number of sites and an approximate translation of it back into numbers:
```
motif: nomptsjiqnnooo
Nsites: 76
[-24.425044699211014, -0.86823247104752, -47.98185692737451, 22.688579757115974, 116.91582866976995, 93.35901644160646, -118.65229361186499, -142.20910584002849, 46.24539198527947, -24.425044699211014, -24.425044699211014, -0.86823247104752, -0.86823247104752, -0.86823247104752]
```

Here are the commands to recreate this for polynomials of degree 1 to 6:

```
./run_part1.sh
./run_part2.sh
./run_part3.sh
```

## Prerequisites
[meme](http://meme-suite.org/doc/download.html)\
Beautiful Soup\
NumPy\
matplotlib\
seaborn\
argparse\
json
