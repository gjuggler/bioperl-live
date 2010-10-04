# $Id: newick.pm 16108 2009-09-16 17:07:49Z cjfields $
#
# BioPerl module for Bio::TreeIO::newick
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::newick - TreeIO implementation for parsing 
  Newick/New Hampshire/PHYLIP format.

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  my $treeio = Bio::TreeIO->new(-format => 'newick', 
                               -file => 't/data/LOAD_Ccd1.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of Newick/PHYLIP/New Hampshire format.

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::TreeIO::newick;
use strict;
use Switch;

use Bio::Event::EventGeneratorI;

use base qw(Bio::TreeIO Bio::TreeIO::NewickParser);

=head2 new

 Title   : new
 Args    : -print_count     => boolean  default is false
           -bootstrap_style => set the bootstrap style (one of nobranchlength,
							molphy, traditional)
           -order_by        => set the order by sort method 
                               (see L<Bio::Node::Node::each_Descendent()> )

=cut

sub _initialize {
  my $self = shift;
  $self->SUPER::_initialize(@_);
  my ($print_count) = $self->_rearrange( [ qw(PRINT_COUNT) ], @_ );

  if ($print_count) {
    $self->param('print_tree_count',1);
  }
  return;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : L<Bio::Tree::TreeI>
 Args    : none


=cut

sub next_tree {
  my ($self) = @_;
  local $/ = ";\n";
  return unless $_ = $self->_readline;

  s/[\r\n]//gs;
  my $score;
  my $despace = sub { my $dirty = shift; $dirty =~ s/\s+//gs; return $dirty };
  my $dequote = sub {
    my $dirty = shift;
    $dirty =~ s/^"?\s*(.+?)\s*"?$/$1/;
    return $dirty;
  };
  s/([^"]*)(".+?")([^"]*)/$despace->($1) . $dequote->($2) . $despace->($3)/egsx;

  if (s/^\s*\[([^\]]+)\]//) {
    my $match = $1;
    $match =~ s/\s//g;
    $match =~ s/lh\=//;
    if ( $match =~ /([-\d\.+]+)/ ) {
      $score = $1;
    }
  }

  $self->_eventHandler->start_document;

  # Call the parse_newick method as defined in NewickParser.pm
  $self->parse_newick($_);

  my $tree = $self->_eventHandler->end_document;

  # Add the tree score afterwards if it exists.
  if ( defined $tree ) {
    $tree->score($score);
    return $tree;
  }
}

=head2 get_default_params

 Title   : get_default_params
 Usage   : my $params = $treeio->get_default_params;
 Function: Returns a hashref containing the default set of parsing & writing parameters for the Newick format.
 Returns : Hashref
 Args    : none

=cut

sub get_default_params {
  my $self = shift;

  return {
    order_by => ''
    , # Changes the tree's node ordering. Only the default value (blank) is currently supported, which leads to alphanumeric sorting of sibling nodes.
    bootstrap_style => 'traditional',    # Can be 'traditional', 'molphy', 'nobranchlength'

    #           'nobranchlength'   --> don't draw any branch lengths - this
    #                                  is helpful if you don't want to have to
    #                                  go through and delete branch len on all nodes
    #           'molphy' --> draw bootstraps (100) like
    #                                  (A:0.11,B:0.22):0.33[100];
    #           'traditional' --> draw bootstraps (100) like
    #                                  (A:0.11,B:0.22)100:0.33;
    internal_node_id => 'id',            # Can be 'id' or 'bootstrap'

    no_branch_lengths       => 0,        # Set to 1 to turn off all output of branch lengths
    no_bootstrap_values     => 0,        # Set to 1 to turn off all output of bootstrap values
    no_internal_node_labels => 0,        # Set to 1 to turn off all output of internal node labels

    newline_each_node => 0,              # Set to 1 to add a new line after each node.
    print_tree_count  => 0,              # Set to 1 to print out the # of trees at the beginning.
  };
}


=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in newick/phylip format
 Returns : none
 Args    : L<Bio::Tree::TreeI> object

=cut

sub write_tree {
  my ( $self, @trees ) = @_;
  if ( $self->param('print_tree_count') ) {
    $self->_print( sprintf( " %d\n", scalar @trees ) );
  }

  my $params = $self->get_params;

  foreach my $tree (@trees) {
    if (!defined $tree
      || ref($tree) =~ /ARRAY/i
      || !$tree->isa('Bio::Tree::TreeI') ) {
      $self->throw("Calling write_tree with non Bio::Tree::TreeI object\n");
    }
    my @data = $self->_write_tree_Helper( $tree->get_root_node, $params );
    $self->_print( join( ',', @data ) . ";" );
  }

  $self->flush if $self->_flush_on_write && defined $self->_fh;
  return;
}

sub _write_tree_Helper {
  my $self = shift;
  my ( $node, $params ) = @_;
  my @data;

  foreach my $n ( $node->each_Descendent( $params->{order_by} ) ) {
    push @data, $self->_write_tree_Helper( $n, $params );
  }

  my $label = $self->_node_as_string( $node, $params );

  if ( scalar(@data) >= 1 ) {
    $data[0] = "(" . $data[0];
    $data[-1] .= ")";
    $data[-1] .= $label;
  } else {
    push @data, $label;
  }

  return @data;
}

=head2 node_as_string

 Title   : node_as_string
 Usage   : my $string = $treeio->_node_as_string($node,$params);
 Function: Returns a string representation of the given node (and *just* the given node -- this is not recursive!)
 Returns : string
 Args    : none

=cut

sub _node_as_string {
  my $self   = shift;
  my $node   = shift;
  my $params = shift;

  my $label_stringbuffer = '';

  if ( $params->{no_bootstrap_values} != 1
    && !$node->is_Leaf
    && defined $node->bootstrap
    && $params->{bootstrap_style}  eq 'traditional'
    && $params->{internal_node_id} eq 'bootstrap' ) {

    # If we're an internal node and we're using 'traditional' bootstrap style,
    # we output the bootstrap instead of any label.
    my $bootstrap = $node->bootstrap;
    $label_stringbuffer .= $bootstrap if ( defined $bootstrap );
  } elsif ( $params->{no_internal_node_labels} != 1 ) {

    # Tack on the internal label unless the $params specifically disallow them.
    my $id = $node->id;
    $label_stringbuffer .= $id if ( defined $id );
  }

  if ( $params->{no_branch_lengths} != 1 ) {

    # Tack on the branch length.
    my $blen = $node->branch_length;
    $label_stringbuffer .= ":" . $blen if ( defined $blen );
  }

  if ( $params->{bootstrap_style} eq 'molphy' && $params->{no_bootstrap_values} != 1 ) {

    # Tack on a molphy-style bootstrap values unless the no_bootstrap_values is set.
    my $bootstrap = $node->bootstrap;
    $label_stringbuffer .= "[$bootstrap]" if ( defined $bootstrap );
  }

  if ( $params->{newline_each_node} == 1 ) {

    # Add a newline if we're meant to.
    $label_stringbuffer .= "\n";
  }

  return $label_stringbuffer;
}

1;
