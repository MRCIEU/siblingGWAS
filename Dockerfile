FROM rocker/r-ver:4.5.2

RUN Rscript -e 'install.packages( \
        c( \
            "data.table", \
            "lmtest", \
            "optparse", \
            "sandwich", \
            "R.utils" \
        ), \
        repos="https://www.stats.bris.ac.uk/R/")'

COPY within-family.R /usr/local/bin/within-family
COPY within-family-plots.R /usr/local/bin/within-family-plots

CMD ["within-family", "-h"]
