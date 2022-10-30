awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
