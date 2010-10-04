# POD documentation - main docs before the code

=head1 NAME

Module which implements a newick string parser as a finite state machine which enables it
to parse the full Newick specification.

Taken largely from the Ensembl Compara file with the same name (Bio::EnsEMBL::Compara::Graph::NewickParser),
this module adapts the parser to work with BioPerl's event handler-based parsing scheme.

This module is used by nhx.pm and newick.pm, and is NOT called directly. Instead, both of those
parsing modules extend this module in order to gain access to the main parsing method.

=head1 SYNOPSIS

  # From newick.pm
  use base qw(Bio::TreeIO Bio::TreeIO::NewickParser);

  # in the next_tree method...
  $self->parse_newick($_);

=head1 DESCRIPTION

This module correctly parses the Newick and NHX formats, sending calls to the BioPerl
TreeEventHandler when appropriate in order to build and populate the node objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jessica Severin (EnsEMBL implementation), Greg Jordan (BioPerl adaptation)

=cut

package Bio::TreeIO::NewickParser;

use strict;
use Switch;

sub parse_newick {
  my $self = shift;
  my $newick = shift;

  $newick = $newick . ";" unless ($newick =~ m/;/);

  my $count=1;
  my $state=1;
  my $bracket_level = 0; # Keeps track of the parser's current depth in the tree.

  my $debug = $self->verbose;

  my $token = next_token(\$newick, "(;");
  $self->_start('tree');

  my $leaf_flag = 0;
  while(defined($token)) {    
    switch ($state) {
      case 1 {
        # A new node.
        $self->_start('node');
        $self->debug("   -> [$token]\n");

        if($token eq '(') {
          # Create a new set.
          $self->debug("    create set\n")  if($debug);
          $token = next_token(\$newick, "[(:,)");
          $state = 1;
          $bracket_level++; # We saw a parenthesis so we're going deeper into the tree...
        } else {
          # No deeper clade, so now we're naming a leaf node.
          $state = 2;
          $leaf_flag = 1;
        }
      }
      case 2 {
        # Naming a node.
        if(!($token =~ /[\[\:\,\)\;]/)) { 
          if (!$leaf_flag && $self->param('internal_node_id') eq 'bootstrap') {
            # We've been told to parse internal node IDs as bootstrap values. Do that here.
            $self->_start('bootstrap');
            $self->_chars($token);
            $self->_end('bootstrap');
            $token = '';
          }

          # Use the token as the ID.
          $self->_start('id');
          $self->_chars($token);
          $self->_end('id');

          $self->debug("    naming leaf\n") if ($debug);
          $token = next_token(\$newick, "[:,);"); # Get the next token and continue on.
        }
        $state = 3;
      }
      case 3 {
        # Optional : and distance
        if($token eq ':') {
          # The token between the colon and the next control char is the branch length.
          $token = next_token(\$newick, "[,);"); 

          $self->_start('branch_length');
          $self->_chars($token);
          $self->_end('branch_length');

          $token = next_token(\$newick, ",);"); # Move to the next control char.
        } elsif ($token eq '[') {
          # Here we've got an NHX tag without a previous branch length.
          $token .= next_token(\$newick, ",);");
        }
        $state = 4;
      }
      case 4 {
        # Optional NHX tags
        if($token =~ /\[\&\&NHX/) {
            $self->_start('nhx_tag');

            # Careful: this regexp gets rid of all NHX wrapping in one step!
            $token =~ /\[\&\&NHX\:(\S+)\]/;
            if ($1) {
                # NHX may be empty, presumably at end of file, just before ";"
                my @attributes = split ':', $1;
                foreach my $attribute (@attributes) {
                    $attribute =~ s/\s+//;
                    my($key,$value) = split '=', $attribute;

                    $self->_start('tag_name');
                    $self->_chars($key);
                    $self->_end('tag_name');

		    $self->_start('tag_value');
		    $self->_chars($value);
		    $self->_end('tag_value');
                  }
            }
            $self->_end('nhx_tag');

            $token = next_token(\$newick, ",);"); # Move to the next control char.
        }
        $state = 5;
      }
      case 5 {
        # End node.
        if($token eq ')') {
          $self->_end('node');

          $token = next_token(\$newick, "[:,);"); # Next control char.
          if (defined $token && $token eq '[') {
            # It is possible to have anonymous internal nodes w/ no name
            # and no blength but with NHX tags
            #
            # We use leaf_flag=0 to let the parser know that it's labeling an internal
            # node. This affects how potential bootstrap values are handled in state 2.
            $leaf_flag = 0;
            $state = 2;
          } else {
            $leaf_flag = 0;
            $state = 2;
          }

          # We saw a right-parenthesis, so we're going one level back up in the tree.
          $bracket_level--;
        } elsif($token eq ',') {
          $self->_end('node');

          $token = next_token(\$newick, "[(:,)"); # Can be branch length-less nhx nodes
          $state=1;
        } elsif($token eq ';') {
          # Done with tree
          die("Parse error: unbalanced ()\n") if($bracket_level ne 0);

          $self->_end('node');
          $self->_end('tree');

          $token = next_token(\$newick, "(");
          $state=13;
        } else {
          $self->debug("[$token]]\n");
          die("Parse error: expected ; or ) or ,\n");
        }
      }
      case 13 {
        die("Parse error: nothing expected after ;");
      }
    }
  }

  if ($self->_eventHandler->within_element('tree')) {
    $self->_end('node');
    $self->_end('tree');
  }
}

sub _chars {
  my $self = shift;
  my $chars = shift;

  $self->_eventHandler->characters($chars);
}

sub _start {
  my $self = shift;
  my $name = shift;

  $self->_eventHandler->start_element({Name=>$name});
}

sub _end {
  my $self = shift;
  my $name = shift;

  $self->_eventHandler->end_element({Name=>$name});
}

sub next_token {
  my $string = shift;
  my $delim = shift;
  
  $$string =~ s/^(\s)+//;

  return undef unless(length($$string));
  
  #print("input =>$$string\n");
  #print("delim =>$delim\n");
  my $index=undef;

  my @delims = split(/ */, $delim);
  foreach my $dl (@delims) {
    my $pos = index($$string, $dl);
    if($pos>=0) {
      $index = $pos unless(defined($index));
      $index = $pos if($pos<$index);
    }
  }
  unless(defined($index)) {
    die("couldn't find delimiter $delim\n $$string");
  }

  my $token ='';

  if($index==0) {
    $token = substr($$string,0,1);
    $$string = substr($$string, 1);
  } else {
    $token = substr($$string, 0, $index);
    $$string = substr($$string, $index);
  }

  #print("  token     =>$token\n");
  #print("  outstring =>$$string\n\n");
  
  return $token;
}



1;
