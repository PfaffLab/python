#!/usr/bin/python


def alterDict(d,szKey):
	if szKey in d:
		d[szKey] += 1
	else:
		d[szKey] = 0


def main():
	dl = {"one":0, "two":0 }

	print "start", dl
	alterDict(dl,"one")
	print "first", dl
	alterDict(dl,"three")
	print "second", dl


main()
