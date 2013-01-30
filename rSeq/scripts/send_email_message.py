import argparse
import base64
import sys

from rSeq.utils.misc import email_notification

def main():
    
    desc  = """ """
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('sender', type=str,
                        help="""what email is this from?""")
    parser.add_argument('to', type=str,
                        help="""what email is this going to?""")
    parser.add_argument('subject', type=str,
                        help="""quoted string. Example: 'important notification'""")
    parser.add_argument('txt', type=str,
                            help="""quoted string. Example: 'this is the body of
                            the important notification i talked about in the subject'""")    
    parser.add_argument('pw_file', type=str,
                            help="""path to file with base64 encoded password.""")    

    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    args = parser.parse_args()

        
    
    pw = open(args.pw_file,'rU').readlines()
    pw = ''.join(pw)
    pw = pw.replace('\n','')
    pw = base64.b64decode(pw)
    email_notification(sender=args.sender,to=args.to,subject=args.subject,txt=args.txt,pw=pw)
    
    
if __name__ == '__main__':
    main()