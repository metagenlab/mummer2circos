#! /usr/bin/python

import re

class Circos_config:
  def __init__(self, caryotype_file,
               chr_spacing_list=[],
               show_ticks="yes",
               show_tick_labels="yes",
               ideogram_spacing=0,
               radius=0.75,
               label_radius=0.05,
               show_ideogram_labels="no",
               color_files=""):

    self.plots = ""
    self.links = ""
    self.highlights = ""
    
    self.template_caryotype= "karyotype = %s\n" \
                             " chromosomes_units           = 10000\n" \
                             " chromosomes_display_default = yes\n" % caryotype_file


    self.template_ideograms = "<ideogram>\n" \
                              " <spacing>\n" \
                              " default            = %su\n" \
                              " %s" \
                              " </spacing>\n" \
                              " \n" \
                              " # thickness and color of ideograms\n" \
                              " thickness          = 12p\n" \
                              " stroke_thickness   = 1\n" \
                              " stroke_color       = black\n" \
                              " \n" \
                              " # the default chromosome color is set here and any value\n" \
                              " # defined in the karyotype file overrides it\n" \
                              " fill               = yes\n" \
                              " fill_color         = black\n" \
                              " \n" \
                              " # fractional radius position of chromosome ideogram within image\n" \
                              " radius             = %sr\n" \
                              " show_label         = %s\n" \
                              " label_font         = default\n" \
                              " label_radius       = dims(ideogram,radius) + %sr\n" \
                              " label_size         = 30\n" \
                              " label_parallel     = no\n" \
                              " \n" \
                              " # show_bands determines whether the outline of cytogenetic bands\n" \
                              " # will be seen\n" \
                              " show_bands         = yes\n" \
                              " band_stroke_thickness = 1\n" \
                              " " \
                              " # in order to fill the bands with the color defined in the karyotype\n" \
                              " # file you must set fill_bands\n" \
                              " fill_bands         = yes\n" \
                              " band_transparency  = 1\n" \
                              " \n" \
                              " </ideogram>\n" % (ideogram_spacing,
                                                  self.add_spacing(chr_spacing_list),
                                                  radius,
                                                  show_ideogram_labels,
                                                  label_radius)

    self.end_ticks = "  <tick>\n" \
                      "  multiplier   = 1\n" \
                      "  position = end\n" \
                      " show_ticks         = yes\n" \
                      "  size              = 4p\n" \
                      "  show_label        = no\n" \
                      "  label_size        = 0p\n" \
                      "  format    = %s bp\n" \
                      "  </tick>\n" \

    self.big_ticks = " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 10u\n" \
                          " size              = 15p\n" \
                          " show_label        = yes\n" \
                          " label_size        = 25p\n" \
                          " format            = %s kb\n" \
                          " thickness         = 2p\n" \
                          " </tick>\n" \


    self.template_ticks = "show_ticks         = %s\n" \
                          " show_tick_labels   = %s\n" \
                          " \n" \
                          " <ticks>\n" \
                          " tick_label_font    = condensed\n" \
                          " radius             = dims(ideogram,radius_outer)\n" \
                          " label_offset       = 8p\n" \
                          " label_size         = 4p\n" \
                          " color              = black\n" \
                          " thickness          = 2p\n" \
                          " \n" \
                          " \n" \
                          " <tick>\n" \
                          " show_ticks         = yes\n" \
                          " skip_first_label = no\n" \
                          " multiplier   = 10/1u\n" \
                          " spacing           = 1u\n" \
                          " size              = 5p\n" \
                          " show_label        = no\n" \
                          " label_size        = 5p\n" \
                          " format            = %s\n" \
                          " </tick>\n" \
                          " %s\n" \
                          " \n" \
                          " \n" \
                          " </ticks>\n" % (show_ticks, show_tick_labels, "%%.1d", self.big_ticks % "%%.1d")

    self.template_rules = "<rules>\n" \
                          " %s\n" \
                          "</rules>\n"

    self.template_backgrounds = "<backgrounds>\n" \
                          " %s\n" \
                          "</backgrounds>\n"

    # #" <<include colors.rn.conf>>\n" \
    self.settings ="<colors>\n" \
                   " %s\n" \
                   " #<<include brewer.all.conf>>\n" \
                   " </colors>\n" \
                   " <image>\n"\
                   " image_map_use      = yes\n" \
                   " image_map_overlay  = no\n" \
                   " image_map_overlay_stroke_color     = red\n" \
                   " <<include image.conf>>\n" \
                   " </image>\n" \
                   " #includes  etc/colors.conf\n" \
                   " #          etc/fonts.conf\n" \
                   " #          etc/patterns.conf\n" \
                   " <<include colors_fonts_patterns.conf>>\n" \
                   " # system and debug settings\n" \
                   " <<include housekeeping.conf>>\n" \
                   " anti_aliasing*     = no\n"

    self.complete_file = self.template_caryotype + self.template_ideograms + self.template_ticks + "%s %s %s" + self.settings % (color_files)

  def _template_spacing(self, chr1, chr2):
    template = '<pairwise %s %s>\n' \
               ' spacing = 2u\n' \
               '</pairwise>\n' % (chr1, chr2)
    return template
    
  def _template_plot(self, file, type="line", r0=1,
                     r1=1.05, color=False, fill_color="red", thickness = "0.8p", z = 1, rules ="", backgrounds="", url="", min=False, max=False):

        template1 = "<plot>\n" \
               "type		    = %s\n" \
               " url                = %s[id]\n" \
               " r0                 = %s\n" \
               " r1                 = %s\n" \
               " fill_color         = %s\n" \
               " thickness          = %s\n" \
               " file               = %s\n" \
               " z                  = %s\n" % (type, url, r0, r1, fill_color, thickness, file, z)

        if color:
            template1+= " color          = %s\n" % color
        if min:
            template1+= " min          = %s\n" % min
        if max:
            max+= " max          = %s\n" % max
        template_rules = " %s\n" \
               " %s\n" \
               " </plot>\n" % (rules, backgrounds)

        return template1 + template_rules

  def template_rule(self, condition, fill_color):
    template = "<rule>\n" \
               " condition          = %s\n" \
               " fill_color         = %s\n" \
               " </rule>\n" % (condition, fill_color)
    return template

  def template_background(self, color):
    template = "<background>\n" \
               " color         = %s\n" \
               " </background>\n" % (color)
    return template


  def _template_link(self, link_file, color="black_a5", thickness=1):
    template = "<link>\n" \
               "ribbon = yes\n" \
               "file          = %s\n" \
               "color         = %s\n" \
               "radius        = 0.99r\n" \
               "bezier_radius = 0.1r\n" \
               "thickness     = %s\n" \
               " </link>\n" % (link_file, color, thickness)
    return template


  def _template_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", url="/chlamdb/locusx/chlamydia_03_15/"):
    template ="<highlight>\n" \
              " fill_color = %s\n" \
              " file       = %s\n" \
              " r1         = %s\n" \
              " r0         = %s\n" \
              " url = %s[id]\n" \
              " </highlight>\n" % (fill_color, file, r1, r0, url)

    return template


  def add_plot(self, file, type="line", r0=1, r1=1.05,
               fill_color="grey_a1", thickness = "2p", z = 1, rules ="", backgrounds="", url="", min=False, max=False, color=False):
    plot = self._template_plot(file, type, r0, r1, color, fill_color, thickness, z, rules, backgrounds, url, min=min, max=max)
    if len(re.findall("</plots>", self.plots))>0:
      # remove end balise
      self.plots = re.sub("</plots>", "", self.plots)
      # add new plot and end balise
      self.plots = self.plots + plot + "</plots>\n"
    else:
      self.plots = "<plots>\n" + plot + "</plots>\n"



  def add_link(self, link_file, color="black_a5", thickness=1):
    link = self._template_link(link_file, color=color, thickness=thickness)
    if len(re.findall("</links>", self.plots))>0:
      # remove end balise
      self.links = re.sub("</links>", "", self.plots)
      # add new plot and end balise
      self.links = self.links + link + "</links>\n"
    else:
      self.links = "<links>\n" + link + "</links>\n"


  def add_highlight(self, file, fill_color="grey_a1", r1="1.55r", r0="1.50r", href=""):
    highlight = self._template_highlight(file, fill_color, r1, r0, href)
    if len(re.findall("</highlights>", self.highlights))>0:
      # remove end balise
      self.highlights = re.sub("</highlights>", "", self.highlights)
      # add new plot and end balise
      self.highlights = self.highlights + highlight + "</highlights>\n"
    else:
      self.highlights = "<highlights>\n" + highlight + "</highlights>\n"

  def add_spacing(self, chr_pair_list):
      all_spacing = ''
      if len(chr_pair_list) == 0:
          return ''
      else:
          for pair in set([tuple(i) for i in chr_pair_list]):
              all_spacing += self._template_spacing(pair[0], pair[1])
          return all_spacing
    
  def get_file(self):
    return self.complete_file % (self.plots, self.highlights, self.links)

