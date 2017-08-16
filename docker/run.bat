docker run 	--rm -it 	^
		-v %CD%\work:/home/pytfa/work 	^
		-v %CD%/..:/src/pytfa		^
		pytfa_docker %*
