for num_cv in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
    do
	cp submit.sub submit_exp2-$num_cv.sub
	echo Rscript bcf/test.r 50 2000 200 2 1 2 $num_cv 1  >> submit_exp2-$num_cv.sub
	qsub submit_exp2-$num_cv.sub
done
