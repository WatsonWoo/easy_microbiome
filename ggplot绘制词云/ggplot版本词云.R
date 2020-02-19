---
title: "ggplot版本词云"
author: "wentao"
date: "2020/2/7"
output: html_document
---

```{css}
 pre code,pre,code {
 white-space:pre!important;
 overflow-x: scroll!important; 
} 
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,
                      fig.width = 7,
                      fig.height = 5,
                      fig.align = "center",
                      warning = FALSE,
                      message = FALSE
                      
                      )
```


# ggplot版本的词云:ggwordcloud


安装R包
```{R}

# 安装
# devtools::install_github("lepennec/ggwordcloud")

```

示例数据

```{R}
data("love_words_small")

head(love_words_small)


set.seed(42)
ggplot(love_words_small, aes(label = word, size = speakers)) +
  geom_text_wordcloud() +
  scale_size_area(max_size = 40) +
  theme_minimal()
```


```{R}

data("love_words")
head(love_words)


set.seed(42)
ggplot( love_words,aes(label = word, size = speakers,color = speakers)) +
  geom_text_wordcloud_area(aes(angle = 45 * sample(-2:2, nrow(love_words),
     replace = TRUE,prob = c(1, 1, 4, 1, 1))),
  mask = png::readPNG(system.file("extdata/hearth.png",
                                  package = "ggwordcloud", mustWork = TRUE
  )),
  rm_outside = TRUE
  ) +
  scale_size_area(max_size = 40) +
  theme_minimal() +
  scale_color_gradient(low = "darkred", high = "red")
#> Some words could not fit on page. They have been removed.


```


```{R}

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE)
set.seed(42)
ggplot(
  love_words_small %>%
    gather(key = "type", value = "speakers", -lang, -word) %>%
    arrange(desc(speakers)),
  aes(label = word, size = speakers)
) +
  geom_text_wordcloud_area() +
  scale_size_area(max_size = 40) +
  theme_minimal() +
  facet_wrap(~type)
```


```{R}

set.seed(42)
ggplot(love_words_small, aes(label = word, size = speakers,
                             label_content = sprintf("%s<span style='font-size:7.5pt'>(%g)</span>", word, speakers))) +
  geom_text_wordcloud_area() +
  scale_size_area(max_size = 40) +
  theme_minimal()

```


```{R}

```


```{R}

```


```{R}

```


```{R}

```


```{R}

```


```{R}

```















