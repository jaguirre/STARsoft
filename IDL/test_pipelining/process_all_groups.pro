pro process_all_groups,groupname,path=path,cd=cd

if ~keyword_set(path) then path = '/scr/starfire/testdata/'
if ~keyword_set(cd) then cd = 'CD004'

fullpath = path+cd+'/'+groupname+'/'

files = file_search(fullpath+groupname+'*.txt')
nfiles = n_elements(files)

for i=0,nfiles-1 do begin
   print, files[i]
   process_one_group,files[i]
endfor

end
