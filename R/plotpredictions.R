plotpredictions<-function(dat=LKsurprisal,
                          maintitle="Predictions of \n the surprisal account",ylabel="reading time"){
  ggplot(data=dat, 
         aes(x=cond, y=pred,group=1)) +
    geom_line(size=1.5) +
    geom_point(size=4) +
    ylab(ylabel) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_discrete(name = NULL,                
                     labels = c("(a)", "(b)", "(c)", "(d)")) +
    scale_y_continuous(name = "reading time",
                       limits = c(450, 750)) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=14, color="black"),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.ticks.x = element_blank()) +
    ggtitle(maintitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", colour = "black")) + magnifytext()
}



