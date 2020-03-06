CoronaVirus HW take 2
================
Rachael McVicar
3/4/2020

## Coronavirus

Here we analyze infection data for the 2019 novel Coronavirus COVID-19
(2019-nCoV) epidemic. The raw data is pulled from the Johns Hopkins
University Center for Systems Science and Engineering (JHU CCSE)
Coronavirus repository.

A CSV file is available here
<https://github.com/RamiKrispin/coronavirus-csv>

``` r
url <- "https://tinyurl.com/COVID-2019"
virus <- read.csv(url)

tail(virus)
```

    ##      Province.State Country.Region      Lat     Long       date cases      type
    ## 2881         Shanxi Mainland China  37.5777 112.2922 2020-03-05     2 recovered
    ## 2882        Sichuan Mainland China  30.6171 102.7103 2020-03-05    19 recovered
    ## 2883        Tianjin Mainland China  39.3054 117.3230 2020-03-05     4 recovered
    ## 2884       Victoria      Australia -37.8136 144.9631 2020-03-05     3 recovered
    ## 2885       Xinjiang Mainland China  41.1129  85.2401 2020-03-05     1 recovered
    ## 2886       Zhejiang Mainland China  29.1832 120.0934 2020-03-05    10 recovered

6 rows and 7 columns

> Q1. How mant total infected cases are there around the world?

``` r
total_cases<- sum(virus$cases)
total_cases
```

    ## [1] 155031

155,031 cases around the world as of 3/6/2020 11:30pm

> Q2. How many deaths linked to infected cases have there been?

Let’s have a look at the *$type* column

``` r
inds<- virus$type == "death"
death_cases<- sum(virus[inds, "cases"])
death_cases
```

    ## [1] 3348

3,348 deaths linked to total cases as of 3/6/2020 11:30pm

> Q3. What is the overall death rate?

percent death is…

``` r
round(death_cases/total_cases * 100, 2)
```

    ## [1] 2.16

2.16% of total COVID19 cases result in death

> Q4. What is the death rate in “Mainland
China”?

``` r
table.China.cases<- subset(virus, Country.Region == "Mainland China", select = TRUE)
total.China.cases<- sum(table.China.cases$cases)
```

cases:
135,675

``` r
table.China.deaths<- subset(table.China.cases, type=="death", select= TRUE) 
total.China.deaths<- sum(table.China.deaths$cases)
```

deaths: 3,013

and for the moment of truth\!

``` r
round(total.China.deaths/total.China.cases * 100, 2)
```

    ## [1] 2.22

2.22% death rate in China. It’s horrible news, but figured the code
out\!

> Q5. What is the death rate in Italy, Iran, and the US? 5A)
ITALY

``` r
table.Italy.cases<- subset(virus, Country.Region == "Italy", select = TRUE)
sum(table.Italy.cases$cases)
```

    ## [1] 4420

``` r
total.Italy.cases<- sum(table.Italy.cases$cases)
```

cases:
4,420

``` r
table.Italy.deaths<- subset(table.Italy.cases, type=="death", select= TRUE) 
sum(table.Italy.deaths$cases)
```

    ## [1] 148

``` r
total.Italy.deaths<- sum(table.Italy.deaths$cases)
```

deaths: 148

``` r
round(total.Italy.deaths/total.Italy.cases * 100, 2)
```

    ## [1] 3.35

3.35% death rate in Italy

> 5B)
IRAN

``` r
table.Iran.cases<- subset(virus, Country.Region == "Iran", select = TRUE)
sum(table.Iran.cases$cases)
```

    ## [1] 4359

``` r
total.Iran.cases<- sum(table.Iran.cases$cases)
```

cases:
4,359

``` r
table.Iran.deaths<- subset(table.Iran.cases, type=="death", select= TRUE) 
sum(table.Iran.deaths$cases)
```

    ## [1] 107

``` r
total.Iran.deaths<- sum(table.Iran.deaths$cases)
```

deaths: 107

``` r
round(total.Iran.deaths/total.Iran.cases * 100, 2)
```

    ## [1] 2.45

2.45% death rate in Iran

> 5C) US

``` r
table.US.cases<- subset(virus, Country.Region == "US", select = TRUE)
sum(table.US.cases$cases)
```

    ## [1] 241

``` r
total.US.cases<- sum(table.US.cases$cases)
```

cases: 241

``` r
table.US.deaths<- subset(table.US.cases, type=="death", select= TRUE) 
sum(table.US.deaths$cases)
```

    ## [1] 12

``` r
total.US.deaths<- sum(table.US.deaths$cases)
```

deaths: 12

``` r
round(total.US.deaths/total.US.cases * 100, 2)
```

    ## [1] 4.98

4.98% death rate in US
