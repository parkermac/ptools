function talk2me(gtag, date_string, backfill, opt_string)
% test of getting arguments from the command line

% write the arguments to a text file
fid = fopen('blah.txt','w');
fprintf(fid,'gtag = %s\n',gtag);
fprintf(fid,'date_string = %s\n',date_string);
fprintf(fid,'backfill = %s\n',backfill);
fprintf(fid,'opt_string = %s\n',opt_string);
fclose(fid);
exit;