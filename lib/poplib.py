#!/usr/bin/python3
############################################################
# For popgen, 10.19
# Common library functions
############################################################

def errorOut(errnum, errmsg, ropt=0):
# Formatting for error messages.
	fullmsg = "| ** Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	if ropt:
		return "\n" + border + "\n" + fullmsg + "\n" + border + "\n";
	else:
		print("\n" + border + "\n" + fullmsg + "\n" + border + "\n");

############################################################

def spacedOut(string, totlen):
# Properly adds spaces to the end of a string to make it a given length
	spaces = " " * (totlen - len(string));
	return string + spaces;

############################################################