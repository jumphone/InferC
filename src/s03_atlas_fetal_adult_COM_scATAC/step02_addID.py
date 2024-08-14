
fo=open('LST_ADULT_ID.txt','w')
i=1
for line in open('LST_ADULT.txt'):
    fo.write('ID_A_'+str(i)+'\t'+line) 
    i=i+1
fo.close()

fo=open('LST_FETAL_ID.txt','w')
i=1
for line in open('LST_FETAL.txt'):
    fo.write('ID_F_'+str(i)+'\t'+line)
    i=i+1
fo.close()


