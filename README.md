# Bitcoin returns time-series analysis
## Paper - **Bitcoin_TV-GARCH_NLP_paper.pdf**

Investors’ attention could be regarded as one of the determinants of Bitcoin dynamic. One of the possible ways to study the impact of investor attention is to extract sentiment measure from public news sources and apply it to the analysis. In light of this, this study presents the analysis of CoinMarketCap Bitcoin data between January 1, 2017 and March 23, 2024, by implementing the TV-GARCH-X model. The extracted sentiment was divided on sentiment of news from traditional sources and news from crypto-specific sources. The results reveal worthy of attention influence of sentiment of news from traditional sources. Precisely, traditional sources sentiment resulted in significant relationship with Bitcoin. Findings suggest that traditional news sources as a source of data for the analysis is more valuable than the crypto-specific one.

## Libraries used: 
1) Basic: numpy, pandas
2) Visualization: mathplotlib
3) NLP: NLTK, spaCy, TextBlob, Gensim
4) Parsing: requests, BeautifulSoup, selenium, fakeagent

## Steps:  
1) Google news parsing using parser in **Google_news_parser.ipynb**
2) Cleaning of data and prepararion for sentiment extraction, first look at data in **Google_news_preprocessing.ipynb**
3) Data preprocessing using NLP techniques (NLTK, spaCy, TextBlob, Gensim) in **Working with text and sentiment (NLP).ipynb**
4) Adding of explanatory variables and final dataset construction in **Dataset construction.ipynb**
5) Statistical analysis and building a models in **TV-GARCH-X.R**
