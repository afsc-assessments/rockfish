
ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: pop nork dusky

pop: pop.tpl 
	@admb -f pop.tpl 

nork: nork.tpl 
	@admb  nork.tpl 

dusky: dusky.tpl 
	@admb  dusky.tpl 

clean:
	@$(RM) *.cpp
	@$(RM) *.htp
	@$(RM) *.obj
