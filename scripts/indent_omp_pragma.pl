#!/usr/bin/env perl

# indent_omp_pragma.pl - indent '#pragma omp' directives to match the
# surrounding code indentation (clang-format cannot do this).
# Other pragmas (#pragma once, #pragma region, ...) are left untouched.
#
# A pragma (e.g. OpenMP) annotates the construct FOLLOWING it, so its
# indentation is taken from the next non-empty, non-preprocessor line.
# Falls back to the previous such line when the pragma is the last
# meaningful line of the file.

use strict;
use warnings;
use Getopt::Long;

my $check = 0;
my $help  = 0;
GetOptions(
    'check' => \$check,
    'help'  => \$help,
) or die "Try '$0 --help'\n";

if ($help || !@ARGV) {
    print <<'USAGE';
Usage: indent_omp_pragma.pl [--check] FILE [FILE...]

Indent C/C++ '#pragma omp' directives to match the following statement's
indentation. Other pragmas are left untouched.
Rewrites files in place (only when changed).

Options:
  --check   do not write; exit 1 if any file would change
  --help    show this message
USAGE
    exit($help ? 0 : 1);
}

# Find the leading whitespace of the nearest reference line, scanning
# from $start in steps of $step (+1 = downwards, -1 = upwards).
sub ref_indent {
    my ($lines, $start, $step) = @_;
    for (my $j = $start; $j >= 0 && $j <= $#$lines; $j += $step) {
        my $ref = $lines->[$j];
        next if $ref =~ /^\s*$/;                                       # empty line
        next if $ref =~ m{^\s*(?://|/\*)};                             # comment line
        next if $ref =~ /^\s*#/ && $ref !~ /^\s*#\s*pragma\s+omp\b/;   # other preprocessor line
        if ($ref =~ /^(\s*)\S/) { return $1; }
    }
    return undef;
}

my $need = 0;
for my $file (@ARGV) {
    open my $fh, '<', $file or die "Cannot open $file: $!";
    my $content = do { local $/; <$fh> };
    close $fh;
    $content = '' unless defined $content;

    my $has_final_nl = ($content eq '' || $content =~ /\n$/);
    my @lines = split /\n/, $content, -1;
    pop @lines if $has_final_nl && @lines;

    my $changed = 0;
    for my $i (0 .. $#lines) {
        my $line = $lines[$i];
        next unless $line =~ /^#\s*pragma\s+omp\b/;
        # align with the next line or previous line
        my $indent = ref_indent(\@lines, $i + 1, 1) // ref_indent(\@lines, $i - 1, -1) // '';
        my $new  = $indent . '#pragma' . substr($line, index($line, 'pragma') + length('pragma'));
        if ($new ne $line) {
            $lines[$i] = $new;
            $changed   = 1;
        }
    }

    if (!$changed) {
        print "Unchanged: $file\n";
        next;
    }
    $need = 1;
    if ($check) {
        print "Would reformat: $file\n";
        next;
    }

    open my $out, '>', $file or die "Cannot write $file: $!";
    my $new_content = join("\n", @lines);
    $new_content .= "\n" if $has_final_nl;
    print $out $new_content;
    close $out;
    print "Formatted: $file\n";
}

exit($check && $need ? 1 : 0);
